/*
 * Lightweight IGV-style tracks rendered on HTML5 canvases.
 *
 * Tracks that share a "group" (same module) share one genomic view
 * (offset + basesPerView), so zooming/panning either track moves both in
 * lockstep. Only the features currently in view are drawn each frame, so
 * interaction stays smooth regardless of sequence length.
 *
 * Two track kinds:
 *   - "seq": the reference DNA sequence (colored bases + a position ruler).
 *   - "mut": mutations placed by genomic position, aligned under the bases.
 *
 * Data is pushed from the Shiny server via the "seqTrackData" and
 * "mutTrackData" custom messages. Clicking a mutation sets a Shiny input
 * (clickInput) so the server can react in the future.
 */
(function () {
  "use strict";

  var BASE_COLORS = {
    A: "#2ECC40", // green
    C: "#0074D9", // blue
    G: "#FF851B", // orange
    T: "#FF4136", // red
    N: "#AAAAAA"  // grey
  };
  var OTHER_COLOR = "#777777";
  // Distinct color per amino acid (stop = black, unknown = grey). Used by the
  // protein track so each residue type is visually identifiable.
  var AA_COLORS = {
    A: "#ff7f00", R: "#4daf4a", N: "#984ea3", D: "#377eb8", C: "#dede00",
    E: "#a65628", Q: "#f781bf", G: "#999999", H: "#66c2a5", I: "#fc8d62",
    L: "#8da0cb", K: "#e78ac3", M: "#e41a1c", F: "#a6d854", P: "#ffd92f",
    S: "#e5c494", T: "#80b1d3", W: "#1b9e77", Y: "#d95f02", V: "#7570b3",
    "*": "#000000", X: "#cccccc"
  };
  // Full amino-acid names (shown on hover in the protein track).
  var AA_NAMES = {
    A: "Alanine", R: "Arginine", N: "Asparagine", D: "Aspartate",
    C: "Cysteine", E: "Glutamate", Q: "Glutamine", G: "Glycine",
    H: "Histidine", I: "Isoleucine", L: "Leucine", K: "Lysine",
    M: "Methionine", F: "Phenylalanine", P: "Proline", S: "Serine",
    T: "Threonine", W: "Tryptophan", Y: "Tyrosine", V: "Valine",
    "*": "Stop", X: "Unknown"
  };
  // Black or white text, whichever reads better on the given background.
  function textColor(hex) {
    var h = hex.replace("#", "");
    var r = parseInt(h.substr(0, 2), 16);
    var gg = parseInt(h.substr(2, 2), 16);
    var b = parseInt(h.substr(4, 2), 16);
    return (0.299 * r + 0.587 * gg + 0.114 * b) > 150 ? "#000" : "#fff";
  }
  var INIT_WINDOW = 60;   // bases visible on first render
  var MIN_BASES = 15;     // most zoomed-in (largest letters)
  var LETTER_MIN_PX = 8;  // draw letters once a base is at least this wide
  var AXIS_H = 18;        // px reserved at the bottom for the position ruler
  var TICK_PX = 90;       // target spacing between axis ticks, in px
  var LEFT_GUTTER = 34;   // px reserved on the left for labels (both tracks)
  // Mutation lanes: A/G/C/T by the changed-to base (both missense and
  // nonsense land on the base they change to), plus "indel" for
  // non-single-base changes.
  var MUT_ROWS = ["A", "G", "C", "T", "indel"];
  var MUT_ROW_INDEX = { A: 0, G: 1, C: 2, T: 3, indel: 4 };

  var groups = {};        // groupId -> shared view + tracks

  function colorFor(b) {
    return BASE_COLORS[b] || OTHER_COLOR;
  }

  // Round a raw spacing up to a "nice" genomic step (1/2/5 x 10^k).
  function niceStep(raw) {
    var pow = Math.pow(10, Math.floor(Math.log10(raw)));
    var f = raw / pow;
    var nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
    return Math.max(1, nf * pow);
  }

  function getGroup(id) {
    var g = groups[id];
    if (!g) {
      g = {
        id: id,
        start: 0, chrom: "", seqLength: 0,
        step: 1,          // genomic coord = start + step * displayIndex
        offset: 0, basesPerView: INIT_WINDOW,
        tracks: {}        // canvasId -> track state
      };
      groups[id] = g;
    }
    return g;
  }

  // Adopt the genomic anchor from an incoming message; reset the view only
  // when the underlying region actually changes (i.e. a new gene).
  function adoptGenome(g, msg) {
    var len = msg.seqLength || (msg.seq ? msg.seq.length : g.seqLength);
    var step = (msg.step === -1) ? -1 : 1;
    var changed = (g.start !== msg.start) || (g.seqLength !== len) ||
                  (g.step !== step);
    g.start = msg.start;
    g.chrom = msg.chrom;
    g.seqLength = len;
    g.step = step;
    if (changed) {
      g.offset = 0;
      g.basesPerView = len || INIT_WINDOW; // default: show the whole gene
    }
  }

  // Display index 0..L-1 reads left->right (gene 5'->3'); genomic coordinate
  // increases (+ strand) or decreases (- strand) with index.
  function genomicAt(g, i) { return g.start + g.step * i; }
  function displayIndexOf(g, v) { return (v - g.start) * g.step; }

  function clampView(g) {
    var n = g.seqLength;
    if (n <= 0) return;
    if (g.basesPerView > n) g.basesPerView = n;
    if (g.basesPerView < MIN_BASES) g.basesPerView = Math.min(MIN_BASES, n);
    if (g.offset < 0) g.offset = 0;
    if (g.offset + g.basesPerView > n) g.offset = Math.max(0, n - g.basesPerView);
  }

  function renderGroup(g) {
    clampView(g);
    for (var id in g.tracks) {
      if (g.tracks.hasOwnProperty(id)) renderTrack(g.tracks[id]);
    }
  }

  // Set device-pixel-ratio-correct size, clear, and return drawing context.
  function prepCanvas(tr) {
    var canvas = tr.canvas;
    var cssW = canvas.clientWidth;
    var cssH = canvas.clientHeight;
    if (!cssW || !cssH) return null; // hidden tab / not laid out yet
    var dpr = window.devicePixelRatio || 1;
    if (canvas.width !== Math.round(cssW * dpr) ||
        canvas.height !== Math.round(cssH * dpr)) {
      canvas.width = Math.round(cssW * dpr);
      canvas.height = Math.round(cssH * dpr);
    }
    var ctx = tr.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    return { cssW: cssW, cssH: cssH, ctx: ctx };
  }

  function drawMessage(p, msg) {
    var ctx = p.ctx;
    ctx.fillStyle = "#888";
    ctx.font = "13px sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(msg || "", p.cssW / 2, p.cssH / 2);
  }

  function renderTrack(tr) {
    var p = prepCanvas(tr);
    if (!p) return;
    var g = tr.group;
    if (tr.empty) { drawMessage(p, tr.message); return; }
    if (g.seqLength <= 0) return;
    if (tr.kind === "seq") drawSeq(tr, p, g);
    else if (tr.kind === "prot") drawProt(tr, p, g);
    else drawMut(tr, p, g);
  }

  // Vertical text, centered in the left gutter (reads bottom-to-top). Font
  // shrinks so the label fits within the available track height (availH).
  function drawVerticalLabel(ctx, text, centerY, availH) {
    var fs = Math.max(7, Math.min(11,
      Math.floor((availH - 10) / (text.length * 0.6))));
    ctx.save();
    ctx.translate(LEFT_GUTTER / 2, centerY);
    ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = "#777";
    ctx.font = fs + "px sans-serif";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(text, 0, 0);
    ctx.restore();
  }

  // Left gutter strip: background, divider line, and a vertical track label.
  function drawGutter(ctx, fillH, label, labelCenterY) {
    ctx.fillStyle = "#f7f7f7";
    ctx.fillRect(0, 0, LEFT_GUTTER, fillH);
    ctx.strokeStyle = "#ccc";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(LEFT_GUTTER + 0.5, 0);
    ctx.lineTo(LEFT_GUTTER + 0.5, fillH);
    ctx.stroke();
    // Available vertical space around the label's center (symmetric).
    var availH = 2 * Math.min(labelCenterY, fillH - labelCenterY);
    drawVerticalLabel(ctx, label, labelCenterY, availH);
  }

  // Protein track: one amino acid per codon (3 bases), aligned to the DNA.
  function isMarkedCodon(tr, k) {
    if (tr.markIndices && tr.markIndices.length) {
      var start = 3 * k;
      var end = start + 3;
      for (var i = 0; i < tr.markIndices.length; i++) {
        var idx = tr.markIndices[i];
        if (idx >= start && idx < end) return true;
      }
      return false;
    }
    return tr.markIndex != null && Math.floor(tr.markIndex / 3) === k;
  }

  function drawProt(tr, p, g) {
    var ctx = p.ctx, cssW = p.cssW, cssH = p.cssH;
    var pw = Math.max(1, cssW - LEFT_GUTTER);
    var bpW = pw / g.basesPerView;
    var codonW = bpW * 3;
    var aa = tr.aa || "";
    var ncod = aa.length;
    var cy = cssH / 2;

    drawGutter(ctx, cssH, "Amino Acid Sequence", cy);

    if (ncod === 0) return;

    var kmin = Math.max(0, Math.floor(g.offset / 3) - 1);
    var kmax = Math.min(ncod - 1, Math.ceil((g.offset + g.basesPerView) / 3));
    var drawLetters = codonW >= LETTER_MIN_PX;

    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, 0, pw, cssH);
    ctx.clip();
    if (drawLetters) {
      var fontPx = Math.min(Math.floor(codonW * 0.7), Math.floor(cssH * 0.6), 16);
      ctx.font = fontPx + "px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
    }
    // Each codon is a cell colored by its amino acid (a protein "barcode" when
    // zoomed out); the residue letter appears once cells are wide enough.
    for (var k = kmin; k <= kmax; k++) {
      var x0 = LEFT_GUTTER + (3 * k - g.offset) * bpW;
      var ch = aa.charAt(k);
      var isMut = isMarkedCodon(tr, k);         // mutated codon: black cell
      var col = isMut ? "#000000" : (AA_COLORS[ch] || AA_COLORS.X);
      ctx.fillStyle = col;
      ctx.fillRect(x0, 0, codonW + 0.5, cssH);
      if (drawLetters) {
        ctx.fillStyle = isMut ? "#ffffff" : textColor(col);
        ctx.fillText(ch, x0 + codonW / 2, cy);
      }
    }
    ctx.restore();
  }

  function isMarkedIndex(tr, i) {
    if (tr.markIndices && tr.markIndices.length) {
      return tr.markIndices.indexOf(i) >= 0;
    }
    return tr.markIndex != null && i === tr.markIndex;
  }

  function drawSeq(tr, p, g) {
    var ctx = p.ctx, cssW = p.cssW, cssH = p.cssH;
    var trackH = Math.max(1, cssH - AXIS_H);   // base blocks occupy the top
    var pw = Math.max(1, cssW - LEFT_GUTTER);  // plot width (right of gutter)
    var bpW = pw / g.basesPerView;             // pixels per base
    var first = Math.floor(g.offset);
    var last = Math.min(tr.seq.length, Math.ceil(g.offset + g.basesPerView));
    var drawLetters = bpW >= LETTER_MIN_PX;

    if (drawLetters) {
      var fontPx = Math.min(Math.floor(bpW * 0.8), Math.floor(trackH * 0.6), 18);
      ctx.font = fontPx + "px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, 0, pw, trackH); // keep bases out of the label gutter
    ctx.clip();
    for (var i = first; i < last; i++) {
      var x = LEFT_GUTTER + (i - g.offset) * bpW;
      var b = tr.seq.charAt(i);
      var isMut = isMarkedIndex(tr, i);
      ctx.fillStyle = isMut ? "#000000" : colorFor(b); // mutated base: black
      ctx.fillRect(x, 0, bpW + 1, trackH); // +1 avoids sub-pixel seams
      if (drawLetters) {
        ctx.fillStyle = "#ffffff";
        ctx.fillText(b, x + bpW / 2, trackH / 2 + 1);
      }
    }
    ctx.restore();

    drawAxis(g, ctx, cssW, trackH, bpW);
    drawGutter(ctx, cssH, "DNA Sequence", trackH / 2);
  }

  // Genomic position ruler; redrawn every frame so it tracks zoom/pan.
  function drawAxis(g, ctx, cssW, trackH, bpW) {
    ctx.fillStyle = "#fafafa";
    ctx.fillRect(0, trackH, cssW, AXIS_H);
    ctx.strokeStyle = "#ccc";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, trackH + 0.5);
    ctx.lineTo(cssW, trackH + 0.5);
    ctx.stroke();

    var tickStep = niceStep((TICK_PX / cssW) * g.basesPerView);
    // Genomic range in view (may run high->low for -1 strand genes).
    var gA = genomicAt(g, g.offset);
    var gB = genomicAt(g, g.offset + g.basesPerView);
    var lo = Math.min(gA, gB), hi = Math.max(gA, gB);
    var firstTick = Math.ceil(lo / tickStep) * tickStep;

    ctx.fillStyle = "#555";
    ctx.strokeStyle = "#999";
    ctx.font = "10px sans-serif";
    ctx.textBaseline = "top";
    ctx.beginPath();
    for (var v = firstTick; v <= hi; v += tickStep) {
      var x = LEFT_GUTTER + (displayIndexOf(g, v) - g.offset) * bpW;
      ctx.moveTo(Math.round(x) + 0.5, trackH);
      ctx.lineTo(Math.round(x) + 0.5, trackH + 4);
      ctx.textAlign = x < 24 ? "left" : (x > cssW - 24 ? "right" : "center");
      ctx.fillText(v.toLocaleString(), x, trackH + 5);
    }
    ctx.stroke();
  }

  function drawMut(tr, p, g) {
    var ctx = p.ctx, cssW = p.cssW, cssH = p.cssH;
    var pw = Math.max(1, cssW - LEFT_GUTTER);
    var bpW = pw / g.basesPerView;
    var laneH = cssH / MUT_ROWS.length;

    tr.hits = [];

    // Horizontal lane separators (across the plot area).
    ctx.strokeStyle = "#eee";
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (var s = 1; s < MUT_ROWS.length; s++) {
      var sy = Math.round(s * laneH) + 0.5;
      ctx.moveTo(LEFT_GUTTER, sy);
      ctx.lineTo(cssW, sy);
    }
    ctx.stroke();

    // Left gutter: A/G/C/T row labels and a divider.
    ctx.fillStyle = "#f7f7f7";
    ctx.fillRect(0, 0, LEFT_GUTTER, cssH);
    ctx.strokeStyle = "#ccc";
    ctx.beginPath();
    ctx.moveTo(LEFT_GUTTER + 0.5, 0);
    ctx.lineTo(LEFT_GUTTER + 0.5, cssH);
    ctx.stroke();
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    for (var li = 0; li < MUT_ROWS.length; li++) {
      var lbl = MUT_ROWS[li];
      ctx.fillStyle = lbl.length === 1 ? colorFor(lbl) : "#777";
      ctx.font = (lbl.length === 1 ? "bold 12px" : "9px") + " monospace";
      ctx.fillText(lbl, LEFT_GUTTER / 2, (li + 0.5) * laneH);
    }

    if (!tr.mutations || tr.mutations.length === 0) {
      ctx.fillStyle = "#aaa";
      ctx.font = "12px sans-serif";
      ctx.fillText("No mutations in this gene", LEFT_GUTTER + pw / 2, cssH / 2);
      return;
    }

    var left = g.offset - 1;
    var right = g.offset + g.basesPerView + 1;
    var r = Math.min(6, laneH / 2 - 2);

    // Markers, clipped to the plot area so they never cover the labels.
    ctx.save();
    ctx.beginPath();
    ctx.rect(LEFT_GUTTER, 0, pw, cssH);
    ctx.clip();
    for (var i = 0; i < tr.mutations.length; i++) {
      var m = tr.mutations[i];
      var ri = MUT_ROW_INDEX[m.base];
      if (ri === undefined) ri = MUT_ROWS.length - 1; // -> "other" lane
      var bi = displayIndexOf(g, m.pos);    // display index along the gene
      if (bi < left || bi > right) continue;
      var x = LEFT_GUTTER + (bi - g.offset + 0.5) * bpW;
      var cy = (ri + 0.5) * laneH;
      var selected = tr.selectedSet && tr.selectedSet[selKey(m)];

      ctx.globalAlpha = selected ? 0.95 : 0.4;
      ctx.fillStyle = selected ? "#000000" : (m.color || OTHER_COLOR);
      ctx.beginPath();
      ctx.arc(x, cy, r, 0, 2 * Math.PI);
      ctx.fill();
      ctx.globalAlpha = 1;

      tr.hits.push({ x: x, y: cy, r: r + 3, m: m });
    }
    ctx.restore();
  }

  function baseAtClientX(tr, clientX) {
    var g = tr.group;
    var rect = tr.canvas.getBoundingClientRect();
    var pw = Math.max(1, rect.width - LEFT_GUTTER);
    var px = clientX - rect.left - LEFT_GUTTER;
    return g.offset + (px / pw) * g.basesPerView;
  }

  function hitAt(tr, clientX, clientY) {
    if (!tr.hits) return null;
    var rect = tr.canvas.getBoundingClientRect();
    var mx = clientX - rect.left, my = clientY - rect.top;
    var best = null, bestD = Infinity;
    for (var i = 0; i < tr.hits.length; i++) {
      var h = tr.hits[i];
      var d = Math.hypot(mx - h.x, my - h.y);
      if (d <= h.r && d < bestD) { best = h; bestD = d; }
    }
    return best;
  }

  function showTip(tr, clientX, clientY, text, multiline) {
    var tip = tr.tip;
    tip.style.whiteSpace = multiline ? "pre" : "nowrap";
    tip.textContent = text;
    tip.style.display = "block";
    // Measure after the content is set, then keep it inside the viewport
    // (flip left / drop below the cursor near an edge).
    var wrapRect = tr.wrap.getBoundingClientRect();
    var tw = tip.offsetWidth, th = tip.offsetHeight;
    var vx = clientX + 12;
    var vy = clientY - th - 6;
    if (vx + tw > window.innerWidth) vx = clientX - tw - 12;
    if (vx < 0) vx = 0;
    if (vy < 0) vy = clientY + 16;
    tip.style.left = (vx - wrapRect.left) + "px";
    tip.style.top = (vy - wrapRect.top) + "px";
  }

  // A mutation is identified by its genomic position + ALT allele (one
  // position can carry several alleles, in different lanes).
  function selKey(m) { return m.pos + "|" + (m.alt || ""); }

  // Single selection: clicking a mutation makes it the (only) selected one;
  // clicking the currently-selected one clears it.
  function toggleSelect(tr, m) {
    var k = selKey(m);
    if (tr.selectedSet && tr.selectedSet[k]) {
      tr.selectedSet = {};
      tr.selectedList = [];
    } else {
      tr.selectedSet = {};
      tr.selectedSet[k] = true;
      tr.selectedList = [m];
    }
  }

  // Send the current (single) selection to the Shiny input, or null if none.
  function pushSelection(tr) {
    if (!tr.clickInput || !window.Shiny) return;
    var m = (tr.selectedList && tr.selectedList[0]) || null;
    Shiny.setInputValue(tr.clickInput, m ? {
      pos: m.pos,
      ref: m.ref,
      alt: m.alt,
      annotation: m.annotation
    } : null, { priority: "event" });
  }

  function escapeHtml(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
      .replace(/>/g, "&gt;");
  }

  // Fill the "selected mutations" box: one line per click, showing the base
  // change, location and annotation. Hidden while the selection is empty.
  function renderSelected(tr) {
    var list = tr.selectedList || [];

    // The action button only appears while a mutation is selected.
    var btn = tr.buttonId && document.getElementById(tr.buttonId);
    if (btn) btn.style.display = list.length ? "inline-block" : "none";

    var box = tr.boxId && document.getElementById(tr.boxId);
    if (!box) return;
    if (list.length === 0) {
      box.style.display = "none";
      box.innerHTML = "";
      return;
    }
    var g = tr.group;
    var html = "<div style='font-weight:bold;margin-bottom:5px;'>" +
      "Selected mutation</div>";
    for (var i = 0; i < list.length; i++) {
      var m = list[i];
      html += "<div>" + escapeHtml(
        g.chrom + ":" + m.pos.toLocaleString() + "  " +
        (m.change || "") + "  (" + (m.annotation || "") + ")"
      ) + "</div>";
    }
    box.innerHTML = html;
    box.style.display = "block";
  }

  function attachHandlers(tr) {
    var canvas = tr.canvas;

    canvas.addEventListener("wheel", function (e) {
      e.preventDefault();
      var g = tr.group;
      if (g.seqLength <= 0) return;
      var anchor = baseAtClientX(tr, e.clientX);
      var factor = e.deltaY < 0 ? 0.85 : 1.18; // up = zoom in
      g.basesPerView *= factor;
      clampView(g);
      var rect = canvas.getBoundingClientRect();
      var pw = Math.max(1, rect.width - LEFT_GUTTER);
      var frac = (e.clientX - rect.left - LEFT_GUTTER) / pw;
      g.offset = anchor - frac * g.basesPerView;
      renderGroup(g);
    }, { passive: false });

    canvas.addEventListener("mousedown", function (e) {
      tr.group.dragging = true;
      tr.group.lastX = e.clientX;
      canvas.style.cursor = "grabbing";
      tr.tip.style.display = "none";
    });

    window.addEventListener("mouseup", function () {
      var g = tr.group;
      if (g.dragging) {
        g.dragging = false;
        canvas.style.cursor = "grab";
      }
    });

    canvas.addEventListener("mousemove", function (e) {
      var g = tr.group;
      if (g.dragging) {
        var rect = canvas.getBoundingClientRect();
        var pw = Math.max(1, rect.width - LEFT_GUTTER);
        var dxBases = ((e.clientX - g.lastX) / pw) * g.basesPerView;
        g.offset -= dxBases;
        g.lastX = e.clientX;
        renderGroup(g);
        return;
      }
      if (tr.empty || g.seqLength <= 0) return;

      if (tr.kind === "seq") {
        var idx = Math.floor(baseAtClientX(tr, e.clientX));
        if (idx < 0 || idx >= tr.seq.length) { tr.tip.style.display = "none"; return; }
        var pos = genomicAt(g, idx);
        showTip(tr, e.clientX, e.clientY,
          g.chrom + ":" + pos.toLocaleString() + "  " + tr.seq.charAt(idx), false);
      } else if (tr.kind === "prot") {
        if (!tr.aa) { tr.tip.style.display = "none"; return; }
        var k = Math.floor(baseAtClientX(tr, e.clientX) / 3);
        if (k < 0 || k >= tr.aa.length) { tr.tip.style.display = "none"; return; }
        var aaCh = tr.aa.charAt(k);
        showTip(tr, e.clientX, e.clientY,
          "residue " + (k + 1) + ": " + (AA_NAMES[aaCh] || aaCh), false);
      } else {
        var h = hitAt(tr, e.clientX, e.clientY);
        if (!h) { tr.tip.style.display = "none"; canvas.style.cursor = "grab"; return; }
        canvas.style.cursor = h.m.clickable ? "pointer" : "default";
        showTip(tr, e.clientX, e.clientY, h.m.hover || "", true);
      }
    });

    canvas.addEventListener("mouseleave", function () {
      tr.tip.style.display = "none";
    });

    if (tr.kind === "mut") {
      canvas.addEventListener("click", function (e) {
        var h = hitAt(tr, e.clientX, e.clientY);
        if (!h || !h.m.clickable) return; // only missense/nonsense
        toggleSelect(tr, h.m);
        renderTrack(tr);
        renderSelected(tr);
        pushSelection(tr);
        hideMutant(tr.group); // a previously generated mutant is now stale
      });
    }

    if (window.ResizeObserver) {
      new ResizeObserver(function () { renderTrack(tr); }).observe(canvas);
    } else {
      window.addEventListener("resize", function () { renderTrack(tr); });
    }
  }

  // Register (once) or fetch a track within its group.
  function getTrack(g, kind, msg) {
    var tr = g.tracks[msg.canvasId];
    if (tr) return tr;
    var canvas = document.getElementById(msg.canvasId);
    if (!canvas) return null;
    tr = {
      group: g,
      kind: kind,
      canvas: canvas,
      ctx: canvas.getContext("2d"),
      tip: document.getElementById(msg.tipId),
      wrap: canvas.parentElement,
      empty: false
    };
    canvas.style.cursor = "grab";
    g.tracks[msg.canvasId] = tr;
    attachHandlers(tr);
    return tr;
  }

  function handleSeq(msg) {
    var g = getGroup(msg.groupId);
    var tr = getTrack(g, "seq", msg);
    if (!tr) return;
    if (msg.empty) {
      tr.empty = true; tr.message = msg.message; tr.seq = "";
      tr.markIndex = null;
      tr.markIndices = [];
      renderTrack(tr);
      return;
    }
    tr.empty = false;
    tr.seq = msg.seq;
    tr.markIndex = (msg.markIndex == null) ? null : msg.markIndex;
    tr.markIndices = Array.isArray(msg.markIndices) ? msg.markIndices.slice() : [];
    adoptGenome(g, msg);
    renderGroup(g);
  }

  function handleMut(msg) {
    var g = getGroup(msg.groupId);
    var tr = getTrack(g, "mut", msg);
    if (!tr) return;
    tr.clickInput = msg.clickInput;
    tr.boxId = msg.boxId;
    tr.buttonId = msg.buttonId;
    g.mutantWrapId = msg.mutantWrapId;
    g.structCompareId = msg.structCompareId;
    // New gene/isolate -> clear the previous selection and any mutant tracks.
    tr.selectedSet = {};
    tr.selectedList = [];
    renderSelected(tr);
    pushSelection(tr);
    hideMutant(g);
    if (msg.empty) {
      tr.empty = true; tr.message = msg.message; tr.mutations = [];
      renderTrack(tr);
      return;
    }
    tr.empty = false;
    tr.mutations = msg.mutations || [];
    adoptGenome(g, msg);
    renderGroup(g);
  }

  function handleProt(msg) {
    var g = getGroup(msg.groupId);
    var tr = getTrack(g, "prot", msg);
    if (!tr) return;
    if (msg.empty) {
      tr.empty = true; tr.message = msg.message; tr.aa = "";
      tr.markIndex = null;
      tr.markIndices = [];
      renderTrack(tr);
      return;
    }
    tr.empty = false;
    tr.aa = msg.aa || "";
    tr.markIndex = (msg.markIndex == null) ? null : msg.markIndex;
    tr.markIndices = Array.isArray(msg.markIndices) ? msg.markIndices.slice() : [];
    adoptGenome(g, msg);
    renderGroup(g);
  }

  function hideMutant(g) {
    if (g && g.mutantWrapId) {
      var w = document.getElementById(g.mutantWrapId);
      if (w) w.style.display = "none";
    }
    // The structure comparison is only relevant once a protein has been
    // copied; a stale mutant (new gene/selection) means it should go away too.
    if (g && g.structCompareId) {
      var sc = document.getElementById(g.structCompareId);
      if (sc) sc.style.display = "none";
    }
  }

  // "Generate Sequence": render mutant DNA + protein tracks (with a mark at
  // the mutation) in the same view group, then reveal the container.
  function handleMutant(msg) {
    var g = getGroup(msg.groupId);
    var dna = getTrack(g, "seq",
      { canvasId: msg.dnaCanvasId, tipId: msg.dnaTipId });
    var dnaMarks = Array.isArray(msg.dnaMarkIndices) ? msg.dnaMarkIndices.slice() :
      (Array.isArray(msg.markIndices) ? msg.markIndices.slice() : []);
    var protMarks = Array.isArray(msg.protMarkIndices) ? msg.protMarkIndices.slice() :
      (Array.isArray(msg.markIndices) ? msg.markIndices.slice() : []);

    if (dna) {
      dna.empty = false;
      dna.seq = msg.seq;
      dna.markIndex = (msg.markIndex == null) ? null : msg.markIndex;
      dna.markIndices = dnaMarks;
    }
    var prot = getTrack(g, "prot",
      { canvasId: msg.protCanvasId, tipId: msg.protTipId });
    if (prot) {
      prot.empty = false;
      prot.aa = msg.aa;
      prot.markIndex = null;
      prot.markIndices = protMarks;
    }
    g.mutantProt = prot; // latest mutant protein, for the copy button

    // Wire the "copy mutated protein" button once. Copies the protein without
    // the stop codon ("*").
    var btn = msg.copyBtnId && document.getElementById(msg.copyBtnId);
    if (btn && !btn._copyWired) {
      btn._copyWired = true;
      btn.addEventListener("click", function () {
        var seq = ((g.mutantProt && g.mutantProt.aa) || "").replace(/\*/g, "");
        if (navigator.clipboard) navigator.clipboard.writeText(seq);
        var label = btn.textContent;
        btn.textContent = "Copied!";
        setTimeout(function () { btn.textContent = label; }, 1500);
        // Now that the user has a sequence to fold, reveal the structure
        // comparison section (text + file upload + the two Mol* viewers).
        if (g.structCompareId) {
          var sc = document.getElementById(g.structCompareId);
          if (sc) sc.style.display = "block";
        }
      });
    }

    var wrap = document.getElementById(msg.wrapId);
    if (wrap) wrap.style.display = "block";
    adoptGenome(g, msg);
    renderGroup(g);
  }

  if (window.Shiny) {
    Shiny.addCustomMessageHandler("seqTrackData", handleSeq);
    Shiny.addCustomMessageHandler("protTrackData", handleProt);
    Shiny.addCustomMessageHandler("mutantTracks", handleMutant);
    Shiny.addCustomMessageHandler("mutTrackData", handleMut);
  }
})();
