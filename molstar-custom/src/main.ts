import { DefaultPluginSpec, PluginSpec } from 'molstar/lib/mol-plugin/spec';
import { PluginContext } from 'molstar/lib/mol-plugin/context';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { Color } from 'molstar/lib/mol-util/color';
import { StateTransforms } from 'molstar/lib/mol-plugin-state/transforms';
import { StructureElement } from 'molstar/lib/mol-model/structure';
import { atoms } from 'molstar/lib/mol-model/structure/query/queries/generators';
import { StructureProperties } from 'molstar/lib/mol-model/structure';
import { QueryContext } from 'molstar/lib/mol-model/structure/query/context';
import { OrderedSet } from 'molstar/lib/mol-data/int';
import { createStructureRepresentationParams } from 'molstar/lib/mol-plugin-state/helpers/structure-representation-params';
import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder';
import { Structure } from 'molstar/lib/mol-model/structure';
import { StructureSelection } from 'molstar/lib/mol-model/structure';
import { alignAndSuperpose } from 'molstar/lib/mol-model/structure/structure/util/superposition';
import { Camera } from 'molstar/lib/mol-canvas3d/camera';

// Defining Shiny as a global object to allow communication from R Shiny
declare global {
  interface Window {
    Shiny?: {
      addCustomMessageHandler: (type: string, handler: (message: any) => void) => void;
      setInputValue: (name: string, value: any, opts?: { priority?: 'event' | 'default' }) => void;
    };
  }
}

// ---------------------------------------------------------------------------
// Multiple viewers, one per canvas.
//
// Each Shiny module (or each canvas within a module) drives its own Mol*
// viewer. Viewers are keyed by canvas id so that a single module can host more
// than one (e.g. the protein-prediction tab). Every message from R carries the
// target `canvasId`, plus, for `initMolstar`, the DOM ids and the exact Shiny
// input names this viewer should report hover events to.
// ---------------------------------------------------------------------------

interface OverpaintLayer {
  bundle: StructureElement.Bundle;
  color: Color;
  clear: boolean;
}

// Where a structure comes from: AlphaFold DB (by UniProt id) or raw file text
// (an uploaded PDB/mmCIF, e.g. an AlphaFold-predicted mutant).
type StructureSource =
  | { kind: 'alphafold'; uniprotId: string }
  | { kind: 'data'; format: 'pdb' | 'mmcif'; data: string };

interface Viewer {
  plugin: PluginContext;
  overpaintLayers: OverpaintLayer[];
  numInput?: string;          // Shiny input id for hovered residue number
  aaInput?: string;           // Shiny input id for hovered residue amino acid
  refLoci?: StructureElement.Loci; // reference (first-loaded) CA atoms, for superposition
}

const viewers = new Map<string, Viewer>();

// State ref for the reference structure's cartoon, so overpaint coloring
// (pLDDT / domain / motif) can target it via .to(cartoonRef).
const cartoonRef = 'cartoon-representation';

const MySpec: PluginSpec = {
  ...DefaultPluginSpec(),
  config: [
    [PluginConfig.VolumeStreaming.Enabled, false]
  ]
};

function getViewer(msg: { canvasId?: string }): Viewer | undefined {
  const v = msg && msg.canvasId ? viewers.get(msg.canvasId) : undefined;
  if (!v) console.warn('Mol*: no viewer for canvas', msg && msg.canvasId);
  return v;
}

// CA-atom loci of a structure, used as the atom set for superposition.
function caLoci(structure: Structure): StructureElement.Loci {
  const query = atoms({
    atomTest: ctx => StructureProperties.atom.label_atom_id(ctx.element) === 'CA'
  });
  const selection = query(new QueryContext(structure));
  return StructureSelection.toLociWithSourceUnits(selection);
}

// Fetch/parse a structure source into a Mol* data node.
async function downloadSource(plugin: PluginContext, source: StructureSource) {
  if (source.kind === 'alphafold') {
    const data = await plugin.builders.data.download(
      { url: 'https://alphafold.ebi.ac.uk/files/AF-' + source.uniprotId + '-F1-model_v6.pdb' },
      { state: { isGhost: true } }
    );
    return { data, format: 'pdb' as const };
  }
  const data = await plugin.builders.data.rawData(
    { data: source.data }, { state: { isGhost: true } }
  );
  return { data, format: source.format };
}

// Load a structure into a viewer and add a cartoon representation. Returns the
// structure state object (has .ref and .obj?.data) so callers can superpose.
async function loadStructure(v: Viewer, source: StructureSource, colorValue: number) {
  const plugin = v.plugin;
  const { data, format } = await downloadSource(plugin, source);
  const trajectory = await plugin.builders.structure.parseTrajectory(
    data, format === 'mmcif' ? 'mmcif' : 'pdb'
  );
  const model = await plugin.builders.structure.createModel(trajectory);
  const structureData = await plugin.builders.structure.createStructure(
    model, { name: 'model', params: {} }
  );
  return { structureData, colorValue };
}

async function addCartoon(
  plugin: PluginContext, structureData: any, colorValue: number, ref?: string
) {
  const polymer = await plugin.builders.structure.tryCreateComponentStatic(structureData, 'polymer');
  if (!polymer) { console.warn('No polymer component found'); return; }
  await plugin.build()
    .to(polymer)
    .apply(StateTransforms.Representation.StructureRepresentation3D, {
      type: { name: 'cartoon', params: {} },
      colorTheme: { name: 'uniform', params: { value: Color(colorValue) } }
    }, ref ? { ref } : undefined)   // ref lets overpaint (.to(cartoonRef)) find it
    .commit();
}

// Create a viewer for a canvas (once) and wire hover -> Shiny inputs. Returns
// undefined if the canvas/parent DOM or WebGL init is unavailable.
async function createViewer(
  canvasId: string, parentId: string, numInput?: string, aaInput?: string
): Promise<Viewer | undefined> {
  const plugin = new PluginContext(MySpec);
  await plugin.init();

  const canvas = document.getElementById(canvasId) as HTMLCanvasElement | null;
  const parent = document.getElementById(parentId) as HTMLDivElement | null;
  if (!canvas || !parent) {
    console.error('Mol*: missing canvas/parent', canvasId, parentId);
    return undefined;
  }
  if (!(await plugin.initViewer(canvas, parent))) {
    console.error('Mol* viewer failed to initialize.');
    return undefined;
  }

  const v: Viewer = { plugin, overpaintLayers: [], numInput, aaInput };
  viewers.set(canvasId, v);

  // Hover reports the residue to THIS viewer's own Shiny inputs.
  plugin.behaviors.interaction.hover.subscribe(e => {
    const info = getResidueInfo(e.current.loci);
    if (!info) return;
    if (v.numInput) window.Shiny!.setInputValue(v.numInput, info.num, { priority: 'event' });
    if (v.aaInput) window.Shiny!.setInputValue(v.aaInput, info.aa, { priority: 'event' });
  });
  return v;
}

// Initialize (or reuse) the viewer for a canvas and load the reference structure.
async function initMolstar(msg: {
  canvasId: string; parentId: string; numInput?: string; aaInput?: string;
  source: StructureSource;
}) {
  let v = viewers.get(msg.canvasId);
  if (!v) {
    v = await createViewer(msg.canvasId, msg.parentId, msg.numInput, msg.aaInput);
    if (!v) return;
  }

  // Fresh reference structure: drop any camera link, clear, load grey.
  unlinkCanvas(msg.canvasId);
  await v.plugin.clear();
  v.overpaintLayers.length = 0;
  v.refLoci = undefined;

  const { structureData } = await loadStructure(v, msg.source, 0xbebebe);
  await addCartoon(v.plugin, structureData, 0xbebebe, cartoonRef);

  // Remember the reference CA atoms so a later structure can be superposed.
  const refStruct = structureData.obj?.data as Structure | undefined;
  if (refStruct) v.refLoci = caLoci(refStruct);

  v.plugin.managers.camera.reset();
  console.log('initMolstar done for', msg.canvasId);
}

// Load the mutant into its OWN viewer, aligned to the WT's orientation so a
// synced camera shows corresponding parts. The WT is ghost-loaded only to
// compute the superposition transform, then removed (never displayed).
async function initAlignedMutant(msg: {
  canvasId: string; parentId: string; source: StructureSource;
  alignUniprot?: string; colorHex?: string; linkWith?: string; markResidue?: number;
  numInput?: string; aaInput?: string;
}) {
  let v = viewers.get(msg.canvasId);
  if (!v) {
    v = await createViewer(msg.canvasId, msg.parentId, msg.numInput, msg.aaInput);
    if (!v) return;
  }

  unlinkCanvas(msg.canvasId);   // drop stale link before reloading
  await v.plugin.clear();
  v.overpaintLayers.length = 0;

  const hex = (msg.colorHex ?? '#d55e00').replace('#', '');
  const colorValue = parseInt(hex, 16);

  const { structureData } = await loadStructure(v, msg.source, colorValue);
  const mobStruct = structureData.obj?.data as Structure | undefined;

  // Align to the WT frame using a ghost-loaded reference, then drop the ref.
  if (mobStruct && msg.alignUniprot) {
    const refData = await v.plugin.builders.data.download(
      { url: 'https://alphafold.ebi.ac.uk/files/AF-' + msg.alignUniprot + '-F1-model_v6.pdb' },
      { state: { isGhost: true } }
    );
    const refTraj = await v.plugin.builders.structure.parseTrajectory(refData, 'pdb');
    const refModel = await v.plugin.builders.structure.createModel(refTraj);
    const refStructData = await v.plugin.builders.structure.createStructure(
      refModel, { name: 'model', params: {} }
    );
    const refStruct = refStructData.obj?.data as Structure | undefined;
    if (refStruct) {
      const results = alignAndSuperpose([caLoci(refStruct), caLoci(mobStruct)]);
      if (results.length) {
        await v.plugin.build()
          .to(structureData.ref)
          .insert(StateTransforms.Model.TransformStructureConformation, {
            transform: { name: 'matrix', params: { data: results[0].bTransform, transpose: false } }
          })
          .commit();
        console.log('align rmsd:', results[0].rmsd);
      }
    }
    await v.plugin.build().delete(refStructData.ref).commit();
  }

  await addCartoon(v.plugin, structureData, colorValue, cartoonRef);
  v.plugin.managers.camera.reset();

  // Black cartoon strip at the mutation residue, on both this (mutant) viewer
  // and the WT viewer. Done here so both structures are already loaded.
  if (msg.markResidue && msg.markResidue > 0) {
    await highlightDomains(v.plugin, v.overpaintLayers, msg.markResidue, msg.markResidue, '#000000');
    const wt = msg.linkWith ? viewers.get(msg.linkWith) : undefined;
    if (wt) await highlightDomains(wt.plugin, wt.overpaintLayers, msg.markResidue, msg.markResidue, '#000000');
  }

  // Link to the WT viewer now that this viewer is fully set up (avoids a race
  // where a separate linkCameras message arrives before the viewer exists).
  if (msg.linkWith) await linkCameras({ canvasA: msg.linkWith, canvasB: msg.canvasId });
}

// Mirror one viewer's camera to another (both directions) so rotating either
// rotates both. Re-linking replaces any previous link for the same pair.
const cameraLinks = new Map<string, () => void>();

// Drop any camera link involving this canvas (call before re-initializing a
// viewer, so intermediate reload states aren't mirrored to the other viewer).
function unlinkCanvas(canvasId: string) {
  for (const [key, teardown] of Array.from(cameraLinks)) {
    if (key.split('|').indexOf(canvasId) !== -1) {
      teardown();
      cameraLinks.delete(key);
    }
  }
}

async function linkCameras(msg: { canvasA: string; canvasB: string }) {
  const a = viewers.get(msg.canvasA);
  const b = viewers.get(msg.canvasB);
  if (!a || !b) { console.warn('Mol*: linkCameras missing viewer'); return; }

  const key = msg.canvasA + '|' + msg.canvasB;
  const prev = cameraLinks.get(key);
  if (prev) prev();

  await Promise.all([a.plugin.canvas3dInitialized, b.plugin.canvas3dInitialized]);
  const cvA = a.plugin.canvas3d;
  const cvB = b.plugin.canvas3d;
  if (!cvA || !cvB) return;

  // Use didDraw (fires every frame a viewer redraws, incl. interactive
  // zoom/rotate/pan) rather than camera.stateChanged, which does not emit on
  // trackball input. On each redraw, copy the source camera to the other if it
  // differs (the equality check prevents an infinite mirror loop).
  const mirror = (src: typeof cvA, dst: typeof cvB) => src.didDraw.subscribe(() => {
    const snap = src.camera.getSnapshot();
    if (Camera.areSnapshotsEqual(dst.camera.getSnapshot(), snap)) return;
    dst.camera.setState(snap, 0);
    dst.requestDraw();
  });

  const subA = mirror(cvA, cvB);
  const subB = mirror(cvB, cvA);
  cameraLinks.set(key, () => { subA.unsubscribe(); subB.unsubscribe(); });
}

// Empty a viewer: drop any camera link and clear its loaded structure(s) and
// overpaint, leaving a blank canvas. Used e.g. when the gene changes so a
// previously uploaded mutant doesn't linger into the next gene.
async function clearViewer(msg: { canvasId: string }) {
  const v = viewers.get(msg.canvasId);
  if (!v) return;
  unlinkCanvas(msg.canvasId);
  await v.plugin.clear();
  v.overpaintLayers.length = 0;
  v.refLoci = undefined;
}

// ---------------------------------------------------------------------------
// Message handlers — every message routes to a viewer by canvasId.
// ---------------------------------------------------------------------------

window.Shiny?.addCustomMessageHandler('initMolstar', (msg: any) => { initMolstar(msg); });

window.Shiny?.addCustomMessageHandler('clearViewer', (msg: any) => { clearViewer(msg); });

window.Shiny?.addCustomMessageHandler('initAlignedMutant', (msg: any) => { initAlignedMutant(msg); });

window.Shiny?.addCustomMessageHandler('linkCameras', (msg: any) => { linkCameras(msg); });

window.Shiny?.addCustomMessageHandler('highlightDomains', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  highlightDomains(v.plugin, v.overpaintLayers, msg.residueStart, msg.residueEnd, msg.colorHex ?? '#ff0000');
});

window.Shiny?.addCustomMessageHandler('highlightResidueWithSphere', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  highlightResidueWithSphere(v.plugin, [msg.positions], msg.colorHex ?? '#ff0000');
});

window.Shiny?.addCustomMessageHandler('applyPlddtColoring', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  applyPlddtColoring(v.plugin, v.overpaintLayers);
});

window.Shiny?.addCustomMessageHandler('zoomToResidue', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  zoomToResidue(v.plugin, msg.residueNumber, msg.chainId ?? 'A');
});

window.Shiny?.addCustomMessageHandler('clearOverlays', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  clearOverlays(v.plugin, v.overpaintLayers);
});

window.Shiny?.addCustomMessageHandler('resetCamera', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  resetCamera(v.plugin);
});

window.Shiny?.addCustomMessageHandler('clearPaint', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  clearPaint(v.plugin, v.overpaintLayers);
});

window.Shiny?.addCustomMessageHandler('clearSpheres', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  clearSpheres(v.plugin);
});

window.Shiny?.addCustomMessageHandler('takeScreenshot', (msg: any) => {
  const v = getViewer(msg);
  if (!v) return;
  takeScreenshot(v.plugin, msg.canvasId);
});

// ---------------------------------------------------------------------------
// Structure helpers (now take an explicit plugin / overpaint cache).
// ---------------------------------------------------------------------------

export async function highlightResidueWithSphere(
  plugin: PluginContext,
  positions: number[],
  colorHex: string = '#ff0000',
) {
  const flatten = (arr: any): any[] => {
    if (!Array.isArray(arr)) return [arr];
    return arr.reduce((acc: any[], v: any) => {
      if (Array.isArray(v)) return acc.concat(flatten(v));
      return acc.concat(v);
    }, []);
  };

  const rawFlat = flatten(positions);
  const posNums = rawFlat
    .map((p: any) => {
      const n = Number(p);
      return Number.isFinite(n) ? Math.floor(n) : null;
    })
    .filter((n: number | null): n is number => n !== null);

  const uniquePos = Array.from(new Set(posNums)).sort((a, b) => a - b);
  if (uniquePos.length === 0) {
    console.warn('No valid positions after cleaning; aborting highlight.');
    return;
  }

  const hex = colorHex.startsWith('#') ? colorHex.slice(1) : colorHex;
  const colorValue = parseInt(hex, 16);
  if (isNaN(colorValue)) { console.warn('Invalid color:', colorHex); return; }

  const structCell = plugin.managers.structure.hierarchy.current.structures[0];

  const b = plugin.build().to(structCell.cell);
  const group = b.apply(StateTransforms.Misc.CreateGroup, { label: 'mutations' }, { ref: 'mutations' });

  const expression = MS.struct.generator.atomGroups({
    'chain-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_asym_id(), 'A']),
    'residue-test': MS.core.set.has([MS.set(...uniquePos), MS.struct.atomProperty.macromolecular.label_seq_id()]),
    'atom-test': MS.core.rel.eq([MS.struct.atomProperty.macromolecular.label_atom_id(), 'CA'])
  });

  group
    .apply(StateTransforms.Model.StructureSelectionFromExpression, { expression })
    .apply(
      StateTransforms.Representation.StructureRepresentation3D,
      createStructureRepresentationParams(plugin, structCell.cell.obj!.data, {
        type: 'ball-and-stick',
        color: 'uniform',
        colorParams: { value: Color(colorValue) },
        size: 'uniform',
        sizeParams: { value: 10 }
      }),
      { tags: ['mutations-group'] }
    );

  await b.commit();
}

function addColorLayer(
  plugin: PluginContext,
  overpaintLayers: OverpaintLayer[],
  structure: Structure,
  residueIds: number[],
  colorHex: string
) {
  const hex = colorHex.startsWith('#') ? colorHex.slice(1) : colorHex;
  const colorValue = parseInt(hex, 16);

  const query = atoms({
    residueTest: ctx => {
      const seqId = StructureProperties.residue.label_seq_id(ctx.element);
      return residueIds.includes(seqId);
    },
  });

  const selection = query(new QueryContext(structure));
  const bundle = StructureElement.Bundle.fromSelection(selection);

  overpaintLayers.push({ bundle, color: Color(colorValue), clear: false });
}

async function clearOverlays(plugin: PluginContext, overpaintLayers: OverpaintLayer[]) {
  overpaintLayers.length = 0;
  await plugin.build()
    .to(cartoonRef)
    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, { layers: [] })
    .commit();
  plugin.build().delete('mutations').commit();
  console.log('Overpaint layers cleared');
}

async function clearPaint(plugin: PluginContext, overpaintLayers: OverpaintLayer[]) {
  overpaintLayers.length = 0;
  await plugin.build()
    .to(cartoonRef)
    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, { layers: [] })
    .commit();
}

export async function clearSpheres(plugin: PluginContext) {
  plugin.build().delete('mutations').commit();
}

async function highlightDomains(
  plugin: PluginContext,
  overpaintLayers: OverpaintLayer[],
  residueStart: number,
  residueEnd: number,
  colorHex: string = '#ff0000'
) {
  const structure = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!structure) { console.warn('No Structure loaded'); return; }

  const hex = colorHex.startsWith('#') ? colorHex.slice(1) : colorHex;
  const colorValue = parseInt(hex, 16);
  if (isNaN(colorValue)) { console.warn('Invalid color:', colorHex); return; }

  const query = atoms({
    residueTest: ctx => {
      const seqId = StructureProperties.residue.label_seq_id(ctx.element);
      return seqId >= residueStart && seqId <= residueEnd;
    },
  });

  const selection = query(new QueryContext(structure));
  const bundle = StructureElement.Bundle.fromSelection(selection);

  overpaintLayers.push({ bundle, color: Color(colorValue), clear: false });

  await plugin.build()
    .to(cartoonRef)
    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, { layers: overpaintLayers })
    .commit();
}

function getResidueInfo(loci: any): { aa: string, num: number, label: string } | undefined {
  if (!StructureElement.Loci.is(loci) || !loci.elements || loci.elements.length === 0) return;

  const e = loci.elements[0];
  if (!e.unit) return;
  const unit = e.unit;
  const localIndex = OrderedSet.start(e.indices);
  if (localIndex == null) return;
  const atomIndex = unit.elements[localIndex];
  if (atomIndex == null) return;

  const model = unit.model;
  const residueIndex = model.atomicHierarchy.residueAtomSegments.index[atomIndex];
  const compId = model.atomicHierarchy.atoms.label_comp_id.value(atomIndex);
  const seqId = model.atomicHierarchy.residues.label_seq_id.value(residueIndex);

  return { aa: String(compId), num: Number(seqId), label: `${compId} ${seqId}` };
}

export async function zoomToResidue(
  plugin: PluginContext,
  residueNumber: number,
  chainId: string = 'A'
) {
  const structure = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!structure) { console.warn('No structure loaded'); return; }

  const query = atoms({
    chainTest: ctx => StructureProperties.chain.label_asym_id(ctx.element) === chainId,
    residueTest: ctx => StructureProperties.residue.label_seq_id(ctx.element) === residueNumber,
  });

  const selection = query(new QueryContext(structure));
  const bundle = StructureElement.Bundle.fromSelection(selection);
  const loci = StructureElement.Bundle.toLoci(bundle, structure);
  plugin.managers.camera.focusLoci(loci);
  plugin.managers.interactivity.lociHighlights.highlight({ loci });
}

export async function resetCamera(plugin: PluginContext) {
  plugin.managers.camera.reset();
}

export async function applyPlddtColoring(plugin: PluginContext, overpaintLayers: OverpaintLayer[]) {
  const structure = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data;
  if (!structure) { console.warn('No structure loaded'); return; }

  const model = structure.models[0];
  const { atomicConformation, atomicHierarchy } = model;
  const bFactors = atomicConformation.B_iso_or_equiv;
  const residueAtomSegments = atomicHierarchy.residueAtomSegments;
  const residueCount = atomicHierarchy.residues._rowCount;

  const veryHighConf: number[] = [];
  const highConf: number[] = [];
  const lowConf: number[] = [];
  const veryLowConf: number[] = [];

  for (let resIndex = 0; resIndex < residueCount; resIndex++) {
    const atomStartIndex = residueAtomSegments.offsets[resIndex];
    const plddt = bFactors.value(atomStartIndex);
    const seqId = atomicHierarchy.residues.label_seq_id.value(resIndex);
    if (plddt > 90) veryHighConf.push(seqId);
    else if (plddt > 70) highConf.push(seqId);
    else if (plddt > 50) lowConf.push(seqId);
    else veryLowConf.push(seqId);
  }

  if (veryHighConf.length > 0) addColorLayer(plugin, overpaintLayers, structure, veryHighConf, '#0053D6');
  if (highConf.length > 0) addColorLayer(plugin, overpaintLayers, structure, highConf, '#65CBF3');
  if (lowConf.length > 0) addColorLayer(plugin, overpaintLayers, structure, lowConf, '#FFDB13');
  if (veryLowConf.length > 0) addColorLayer(plugin, overpaintLayers, structure, veryLowConf, '#FF7D45');

  await plugin.build()
    .to(cartoonRef)
    .apply(StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle, { layers: overpaintLayers })
    .commit();
}

export async function takeScreenshot(
  plugin: PluginContext,
  canvasId: string,
  resolution: number = 5
) {
  try {
    if (!plugin.canvas3d) { console.error('Canvas3D not available'); return; }
    const canvas = document.getElementById(canvasId) as HTMLCanvasElement;
    const width = canvas.width * resolution;
    const height = canvas.height * resolution;

    const imagePass = plugin.canvas3d.getImagePass({
      multiSample: { mode: 'on', sampleLevel: 2 },
      transparentBackground: false,
      postprocessing: plugin.canvas3d.props.postprocessing
    });

    imagePass.setSize(width, height);
    const imageData = await imagePass.getImageData(plugin.canvas3d.webgl as any, width, height);

    const tempCanvas = document.createElement('canvas');
    tempCanvas.width = width;
    tempCanvas.height = height;
    const ctx = tempCanvas.getContext('2d');
    if (ctx) {
      ctx.putImageData(imageData, 0, 0);
      tempCanvas.toBlob((blob) => {
        if (!blob) return;
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = `molstar_screenshot_${Date.now()}.png`;
        link.click();
        URL.revokeObjectURL(url);
      }, 'image/png');
    }
  } catch (error) {
    console.error('Screenshot failed:', error);
  }
}
