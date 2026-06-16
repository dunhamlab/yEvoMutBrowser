/*
 * Lightweight IGV-style DNA sequence track rendered on an HTML5 canvas.
 *
 * Only the bases currently in view are drawn each frame, so zoom/pan stays
 * smooth regardless of total sequence length. Genomic coordinates are shown
 * in a hover tooltip rather than on an axis.
 *
 * Data is pushed from the Shiny server via the "seqTrackData" custom message.
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
  var INIT_WINDOW = 60;   // bases visible on first render
  var MIN_BASES = 15;     // most zoomed-in (largest letters)
  var LETTER_MIN_PX = 8;  // draw letters once a base is at least this wide
  var AXIS_H = 18;        // px reserved at the bottom for the position ruler
  var TICK_PX = 90;       // target spacing between axis ticks, in px

  // Per-canvas state, keyed by canvas element id.
  var registry = {};

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

  function clampView(st) {
    var n = st.seq.length;
    if (st.basesPerView > n) st.basesPerView = n;
    if (st.basesPerView < MIN_BASES) {
      st.basesPerView = Math.min(MIN_BASES, n);
    }
    if (st.offset < 0) st.offset = 0;
    if (st.offset + st.basesPerView > n) {
      st.offset = Math.max(0, n - st.basesPerView);
    }
  }

  function render(st) {
    var canvas = st.canvas;
    var cssW = canvas.clientWidth;
    var cssH = canvas.clientHeight;
    if (!cssW || !cssH) return; // hidden tab / not laid out yet

    var dpr = window.devicePixelRatio || 1;
    if (canvas.width !== Math.round(cssW * dpr) ||
        canvas.height !== Math.round(cssH * dpr)) {
      canvas.width = Math.round(cssW * dpr);
      canvas.height = Math.round(cssH * dpr);
    }
    var ctx = st.ctx;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);

    if (st.empty) {
      ctx.fillStyle = "#666";
      ctx.font = "13px sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(st.message || "No sequence", cssW / 2, cssH / 2);
      return;
    }

    clampView(st);
    var trackH = Math.max(1, cssH - AXIS_H);   // base blocks occupy the top
    var bpW = cssW / st.basesPerView;          // pixels per base
    var first = Math.floor(st.offset);
    var last = Math.min(st.seq.length, Math.ceil(st.offset + st.basesPerView));
    var drawLetters = bpW >= LETTER_MIN_PX;

    if (drawLetters) {
      var fontPx = Math.min(Math.floor(bpW * 0.8), Math.floor(trackH * 0.6), 18);
      ctx.font = fontPx + "px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
    }

    for (var i = first; i < last; i++) {
      var x = (i - st.offset) * bpW;
      var b = st.seq.charAt(i);
      ctx.fillStyle = colorFor(b);
      // +1 to avoid sub-pixel seams between adjacent blocks
      ctx.fillRect(x, 0, bpW + 1, trackH);
      if (drawLetters) {
        ctx.fillStyle = "#ffffff";
        ctx.fillText(b, x + bpW / 2, trackH / 2 + 1);
      }
    }

    drawAxis(st, ctx, cssW, trackH, bpW);
  }

  // Genomic position ruler drawn along the bottom; redrawn every frame so it
  // tracks zoom/pan. Ticks land on "nice" genomic coordinates.
  function drawAxis(st, ctx, cssW, trackH, bpW) {
    ctx.fillStyle = "#fafafa";
    ctx.fillRect(0, trackH, cssW, AXIS_H);
    ctx.strokeStyle = "#ccc";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, trackH + 0.5);
    ctx.lineTo(cssW, trackH + 0.5);
    ctx.stroke();

    var step = niceStep((TICK_PX / cssW) * st.basesPerView);
    var leftGenome = st.start + st.offset;
    var rightGenome = st.start + st.offset + st.basesPerView;
    var firstTick = Math.ceil(leftGenome / step) * step;

    ctx.fillStyle = "#555";
    ctx.strokeStyle = "#999";
    ctx.font = "10px sans-serif";
    ctx.textBaseline = "top";
    ctx.beginPath();
    for (var g = firstTick; g <= rightGenome; g += step) {
      var x = (g - st.start - st.offset) * bpW;
      ctx.moveTo(Math.round(x) + 0.5, trackH);
      ctx.lineTo(Math.round(x) + 0.5, trackH + 4);
      // keep edge labels from spilling outside the canvas
      ctx.textAlign = x < 24 ? "left" : (x > cssW - 24 ? "right" : "center");
      ctx.fillText(g.toLocaleString(), x, trackH + 5);
    }
    ctx.stroke();
  }

  function baseAtClientX(st, clientX) {
    var rect = st.canvas.getBoundingClientRect();
    var px = clientX - rect.left;
    return st.offset + (px / rect.width) * st.basesPerView;
  }

  function attachHandlers(st) {
    var canvas = st.canvas;

    canvas.addEventListener("wheel", function (e) {
      e.preventDefault();
      var anchor = baseAtClientX(st, e.clientX);
      var factor = e.deltaY < 0 ? 0.85 : 1.18; // up = zoom in
      st.basesPerView *= factor;
      clampView(st);
      // keep the base under the cursor fixed
      var rect = canvas.getBoundingClientRect();
      var frac = (e.clientX - rect.left) / rect.width;
      st.offset = anchor - frac * st.basesPerView;
      clampView(st);
      render(st);
    }, { passive: false });

    canvas.addEventListener("mousedown", function (e) {
      st.dragging = true;
      st.lastX = e.clientX;
      canvas.style.cursor = "grabbing";
      st.tip.style.display = "none";
    });

    window.addEventListener("mouseup", function () {
      if (st.dragging) {
        st.dragging = false;
        canvas.style.cursor = "grab";
      }
    });

    canvas.addEventListener("mousemove", function (e) {
      if (st.dragging) {
        var rect = canvas.getBoundingClientRect();
        var dxBases = ((e.clientX - st.lastX) / rect.width) * st.basesPerView;
        st.offset -= dxBases;
        st.lastX = e.clientX;
        clampView(st);
        render(st);
        return;
      }
      if (st.empty) return;
      var idx = Math.floor(baseAtClientX(st, e.clientX));
      if (idx < 0 || idx >= st.seq.length) {
        st.tip.style.display = "none";
        return;
      }
      var pos = st.start + idx;
      st.tip.textContent = st.chrom + ":" + pos.toLocaleString() +
        "  " + st.seq.charAt(idx);
      var wrapRect = st.wrap.getBoundingClientRect();
      st.tip.style.left = (e.clientX - wrapRect.left + 12) + "px";
      st.tip.style.top = (e.clientY - wrapRect.top - 28) + "px";
      st.tip.style.display = "block";
    });

    canvas.addEventListener("mouseleave", function () {
      st.tip.style.display = "none";
    });

    // Redraw when the canvas gains size (tab shown) or the window resizes.
    if (window.ResizeObserver) {
      new ResizeObserver(function () { render(st); }).observe(canvas);
    } else {
      window.addEventListener("resize", function () { render(st); });
    }
  }

  function getState(canvasId, tipId) {
    var st = registry[canvasId];
    if (st) return st;
    var canvas = document.getElementById(canvasId);
    if (!canvas) return null;
    st = {
      canvas: canvas,
      ctx: canvas.getContext("2d"),
      tip: document.getElementById(tipId),
      wrap: canvas.parentElement,
      seq: "",
      start: 0,
      chrom: "",
      offset: 0,
      basesPerView: INIT_WINDOW,
      dragging: false,
      empty: false
    };
    canvas.style.cursor = "grab";
    registry[canvasId] = st;
    attachHandlers(st);
    return st;
  }

  function handle(msg) {
    var st = getState(msg.canvasId, msg.tipId);
    if (!st) return;
    if (msg.empty) {
      st.empty = true;
      st.message = msg.message;
      st.seq = "";
      render(st);
      return;
    }
    st.empty = false;
    st.seq = msg.seq;
    st.start = msg.start;
    st.chrom = msg.chrom;
    st.offset = 0;
    st.basesPerView = Math.min(INIT_WINDOW, msg.seq.length);
    render(st);
  }

  if (window.Shiny) {
    Shiny.addCustomMessageHandler("seqTrackData", handle);
  }
})();
