const fg = "#20252b",
  grid = "#d5dce1",
  gold = "#b45c00",
  teal = "#007c78";
let exportDensity=0;
function addExport(plot,heatmap) {
  const button=document.createElement('button');button.type='button';
  button.className='plot-export';button.textContent='PNG ↓';button.title='Download figure on white background at 4× resolution';
  plot.root.append(button);
  button.onclick=()=>{
    if(!plot.config)return;
    try {
      exportDensity=4;plot.draw();
      const canvas=document.createElement('canvas');
      canvas.width=Math.round(plot.root.clientWidth*4);canvas.height=Math.round(plot.root.clientHeight*4);
      const context=canvas.getContext('2d');context.fillStyle='#fff';context.fillRect(0,0,canvas.width,canvas.height);
      if(heatmap){context.drawImage(plot.pixels,74*4,37*4);context.drawImage(plot.axes,0,0);}
      else context.drawImage(plot.canvas,0,0);
      canvas.toBlob(blob=>{
        if(!blob)return;
        const url=URL.createObjectURL(blob),link=document.createElement('a');
        link.href=url;link.download=`fmcw-${plot.root.id}.png`;link.click();
        setTimeout(()=>URL.revokeObjectURL(url),1000);
      });
    } finally {exportDensity=0;plot.draw();}
  };
}
export const db = (x) => 10 * Math.log10(Math.max(x, 1e-30));
const stops = [
  [0.0, [9, 12, 27]],
  [0.25, [68, 28, 94]],
  [0.5, [164, 46, 87]],
  [0.75, [237, 115, 78]],
  [1, [249, 237, 164]],
];
function magma(t) {
  t = Math.max(0, Math.min(1, t));
  let k = 1;
  while (stops[k][0] < t && k < 4) k++;
  const a = stops[k - 1],
    b = stops[k],
    w = (t - a[0]) / (b[0] - a[0]);
  return a[1].map((v, i) => v * (1 - w) + b[1][i] * w);
}
function diverging(t) {
  t = Math.max(0, Math.min(1, t));
  const a = t < 0.5 ? [37, 91, 139] : [242, 242, 231],
    b = t < 0.5 ? [242, 242, 231] : [189, 74, 40],
    w = t < 0.5 ? t * 2 : (t - 0.5) * 2;
  return a.map((v, i) => v * (1 - w) + b[i] * w);
}
function sizeCanvas(canvas, width, height) {
  const dpr = exportDensity || Math.min(devicePixelRatio || 1, 2);
  canvas.width = Math.round(width * dpr);
  canvas.height = Math.round(height * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}
export class Heatmap {
  constructor(id, onSelect) {
    this.root = document.getElementById(id);
    this.pixels = document.createElement("canvas");
    this.pixels.className = "pixels";
    this.axes = document.createElement("canvas");
    this.axes.className = "axes";
    this.root.append(this.pixels, this.axes);
    const gl = (this.gl = this.pixels.getContext("webgl2", {
      antialias: false,
      preserveDrawingBuffer: true,
    }));
    if (!gl) throw Error("WebGL2 is required for the data maps.");
    const vs = `#version 300 es\nin vec2 p;out vec2 uv;void main(){uv=(p+1.)/2.;gl_Position=vec4(p,0.,1.);}`;
    const fs = `#version 300 es\nprecision highp float;uniform sampler2D image;uniform float lo,hi;uniform int kind;uniform vec4 sourceRect;in vec2 uv;out vec4 color;vec3 magma(float x){vec3 a=vec3(9,12,27)/255.;vec3 b=vec3(68,28,94)/255.;vec3 c=vec3(164,46,87)/255.;vec3 d=vec3(237,115,78)/255.;vec3 e=vec3(249,237,164)/255.;return x<.25?mix(a,b,x*4.):x<.5?mix(b,c,(x-.25)*4.):x<.75?mix(c,d,(x-.5)*4.):mix(d,e,(x-.75)*4.);}void main(){float v=texture(image,mix(sourceRect.xy,sourceRect.zw,uv)).r;if(isnan(v)||isinf(v)){color=vec4(.24,.29,.32,1.);return;}float x=clamp((v-lo)/(hi-lo),0.,1.);vec3 c=kind==0?(x<.5?mix(vec3(37,91,139)/255.,vec3(242,242,231)/255.,x*2.):mix(vec3(242,242,231)/255.,vec3(189,74,40)/255.,(x-.5)*2.)):magma(x);color=vec4(c,1.);}`;
    const compile = (type, source) => {
      const shader = gl.createShader(type);
      gl.shaderSource(shader, source);
      gl.compileShader(shader);
      if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS))
        throw Error(gl.getShaderInfoLog(shader));
      return shader;
    };
    const program = (this.program = gl.createProgram());
    gl.attachShader(program, compile(gl.VERTEX_SHADER, vs));
    gl.attachShader(program, compile(gl.FRAGMENT_SHADER, fs));
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS))
      throw Error(gl.getProgramInfoLog(program));
    gl.useProgram(program);
    const buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]),
      gl.STATIC_DRAW,
    );
    const location = gl.getAttribLocation(program, "p");
    gl.enableVertexAttribArray(location);
    gl.vertexAttribPointer(location, 2, gl.FLOAT, false, 0, 0);
    this.texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.texture);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    this.pixels.onclick = (e) => {
      if(this.suppressClick){this.suppressClick=false;return;}
      if (!this.config || !onSelect) return;
      const r = this.pixels.getBoundingClientRect();
      onSelect(
        this.config.x0 +
          ((e.clientX - r.left) / r.width) * (this.config.x1 - this.config.x0),
      );
    };
    this.enableZoom();
    addExport(this,true);
    new ResizeObserver(() => this.draw()).observe(this.root);
  }
  enableZoom(options = {}) {
    this.zoomOptions = options;
    if(this.zoomBox) {this.updateResetButton();return;}
    this.zoomBox=document.createElement('div');
    this.zoomBox.className='plot-zoom-box';this.zoomBox.hidden=true;
    this.root.append(this.zoomBox);
    this.resetButton=document.createElement('button');
    this.resetButton.type='button';this.resetButton.className='plot-reset';
    this.resetButton.textContent='Reset plot view';this.resetButton.hidden=true;
    this.resetButton.onclick=()=>this.resetZoom();this.root.after(this.resetButton);
    this.pixels.tabIndex=0;
    this.pixels.setAttribute('aria-label','Drag a region to zoom; click to select a chirp');
    const point=event=>{
      const r=this.pixels.getBoundingClientRect();
      return [Math.max(0,Math.min(1,(event.clientX-r.left)/r.width)),Math.max(0,Math.min(1,(event.clientY-r.top)/r.height))];
    };
    this.pixels.addEventListener('pointerdown',event=>{
      if(event.button!==0 || !this.config || this.zoomOptions.enabled?.()===false)return;
      event.preventDefault();window.getSelection()?.removeAllRanges();
      this.pixels.focus({preventScroll:true});this.suppressClick=false;
      this.brush={pointer:event.pointerId,start:point(event),end:point(event)};
      document.documentElement.classList.add('plot-window-dragging');
      this.pixels.setPointerCapture(event.pointerId);
    });
    this.pixels.addEventListener('pointermove',event=>{
      if(!this.brush || this.brush.pointer!==event.pointerId)return;
      event.preventDefault();this.brush.end=point(event);
      const [x,y]=this.brush.start,[xx,yy]=this.brush.end;
      const r=this.pixels.getBoundingClientRect();
      this.zoomBox.hidden=false;
      Object.assign(this.zoomBox.style,{left:`${74+Math.min(x,xx)*r.width}px`,top:`${37+Math.min(y,yy)*r.height}px`,width:`${Math.abs(x-xx)*r.width}px`,height:`${Math.abs(y-yy)*r.height}px`});
    });
    const finish=(event,cancelled=false)=>{
      const brush=this.brush;
      if(!brush || event.pointerId!==brush.pointer)return;
      this.brush=null;this.zoomBox.hidden=true;
      document.documentElement.classList.remove('plot-window-dragging');
      if(this.pixels.hasPointerCapture(event.pointerId))this.pixels.releasePointerCapture(event.pointerId);
      const rect=this.pixels.getBoundingClientRect();
      const dx=Math.abs(brush.start[0]-brush.end[0])*rect.width,dy=Math.abs(brush.start[1]-brush.end[1])*rect.height;
      if(cancelled){this.suppressClick=true;return;}
      if(dx<6 && dy<6)return;
      this.suppressClick=true;
      const c=this.config;
      const xs=dx<6?[0,1]:[Math.min(brush.start[0],brush.end[0]),Math.max(brush.start[0],brush.end[0])];
      const ys=dy<6?[0,1]:[1-Math.max(brush.start[1],brush.end[1]),1-Math.min(brush.start[1],brush.end[1])];
      const bounds={x0:c.x0+xs[0]*(c.x1-c.x0),x1:c.x0+xs[1]*(c.x1-c.x0),y0:c.y0+ys[0]*(c.y1-c.y0),y1:c.y0+ys[1]*(c.y1-c.y0)};
      if(this.zoomOptions.onZoom)this.zoomOptions.onZoom(bounds);
      else this.zoomTo(bounds);
    };
    this.pixels.addEventListener('pointerup',event=>finish(event));
    this.pixels.addEventListener('pointercancel',event=>finish(event,true));
    this.pixels.addEventListener('lostpointercapture',event=>finish(event,true));
    this.pixels.addEventListener('keydown',event=>{
      if(event.key==='Escape' && this.brush){event.preventDefault();finish({pointerId:this.brush.pointer},true);}
    });
  }
  updateResetButton() {
    if(this.resetButton)this.resetButton.hidden=this.zoomOptions.reset===false || !this.viewport;
  }
  zoomTo(bounds) {
    if(!this.baseConfig)return;
    this.viewport={...this.viewport,...bounds};
    this.config={...this.baseConfig,...this.viewport};
    this.updateResetButton();this.draw();
  }
  resetZoom() {
    this.viewport=null;
    if(this.baseConfig)this.config={...this.baseConfig};
    this.updateResetButton();this.draw();
  }
  enableWindow(options) {
    const band = document.createElement("div");
    (this.windows ??= []).push({ band, options });
    band.className = `background-window ${options.tone ?? ""}`;
    const label = document.createElement("button");
    label.type = "button";
    label.className = "window-label";
    label.textContent = options.label ?? "Background";
    label.setAttribute(
      "aria-label",
      options.moveLabel ?? "Move background window",
    );
    band.append(label);
    for (const edge of options.resizable === false ? [] : ["start", "stop"]) {
      const handle = document.createElement("button");
      handle.type = "button";
      handle.className = `window-handle ${edge}`;
      handle.setAttribute("aria-label", `Resize background ${edge}`);
      handle.dataset.edge = edge;
      band.append(handle);
    }
    this.root.append(band);
    band.addEventListener("dragstart", (event) => event.preventDefault());
    const at = (event) => {
      const rect = this.pixels.getBoundingClientRect();
      return (
        this.config.x0 +
        ((event.clientX - rect.left) / rect.width) *
          (this.config.x1 - this.config.x0)
      );
    };
    band.addEventListener("pointerdown", (event) => {
      if (event.button !== 0 || !options.enabled()) return;
      event.preventDefault();
      event.stopPropagation();
      window.getSelection()?.removeAllRanges();
      document.documentElement.classList.add("plot-window-dragging");
      const [start, stop] = options.getRange();
      this.windowDrag = {
        pointer: event.pointerId,
        anchor: at(event),
        start,
        stop,
        edge: event.target.dataset.edge,
      };
      band.setPointerCapture(event.pointerId);
      band.classList.add("dragging");
    });
    band.addEventListener("pointermove", (event) => {
      const drag = this.windowDrag;
      if (!drag || drag.pointer !== event.pointerId) return;
      event.preventDefault();
      const delta = Math.round(at(event) - drag.anchor);
      let start = drag.start,
        stop = drag.stop;
      if (drag.edge === "start")
        start = Math.max(0, Math.min(stop - 1, start + delta));
      else if (drag.edge === "stop")
        stop = Math.min(options.limit(), Math.max(start + 1, stop + delta));
      else {
        start = Math.max(
          0,
          Math.min(options.limit() - (stop - start), start + delta),
        );
        stop = start + drag.stop - drag.start;
      }
      options.preview(start, stop);
    });
    const finish = (event, cancelled = false) => {
      const drag = this.windowDrag;
      if (!drag || event.pointerId !== drag.pointer) return;
      this.windowDrag = null;
      document.documentElement.classList.remove("plot-window-dragging");
      band.classList.remove("dragging");
      if (band.hasPointerCapture(event.pointerId))
        band.releasePointerCapture(event.pointerId);
      if (cancelled) options.preview(drag.start, drag.stop);
      else {
        const current = options.getRange();
        if (current[0] !== drag.start || current[1] !== drag.stop)
          options.commit();
      }
    };
    band.addEventListener("pointerup", (event) => finish(event));
    band.addEventListener("pointercancel", (event) => finish(event, true));
    band.addEventListener("lostpointercapture", (event) => finish(event, true));
    band.addEventListener("keydown", (event) => {
      if (
        !["ArrowLeft", "ArrowRight"].includes(event.key) ||
        !options.enabled()
      )
        return;
      event.preventDefault();
      const delta =
        (event.key === "ArrowLeft" ? -1 : 1) * (event.shiftKey ? 10 : 1);
      let [start, stop] = options.getRange();
      const edge = event.target.dataset.edge;
      if (edge === "start")
        start = Math.max(0, Math.min(stop - 1, start + delta));
      else if (edge === "stop")
        stop = Math.min(options.limit(), Math.max(start + 1, stop + delta));
      else {
        const width = stop - start;
        start = Math.max(0, Math.min(options.limit() - width, start + delta));
        stop = start + width;
      }
      options.preview(start, stop);
      options.commit();
    });
  }
  drawWindow() {
    if (!this.windows || !this.config) return;
    for (const { band, options } of this.windows) {
      const [start, stop] = options.getRange();
      const c = this.config,
        rect = this.pixels.getBoundingClientRect();
      const a = Math.max(start, c.x0),
        b = Math.min(stop, c.x1);
      band.hidden = b <= a;
      if (b <= a) continue;
      const width = this.root.clientWidth - 176;
      Object.assign(band.style, {
        left: `${74 + ((a - c.x0) / (c.x1 - c.x0)) * width}px`,
        top: "37px",
        width: `${Math.max(options.resizable === false ? 3 : 0, ((b - a) / (c.x1 - c.x0)) * width)}px`,
        height: `${rect.height}px`,
      });
      const label = band.querySelector(".window-label");
      label.textContent = options.getLabel?.() ?? options.label ?? "Background";
      const labelCenter = Math.max(
        label.offsetWidth / 2,
        Math.min(
          width - label.offsetWidth / 2,
          (((a + b) / 2 - c.x0) / (c.x1 - c.x0)) * width,
        ),
      );
      label.style.left = `${labelCenter - ((a - c.x0) / (c.x1 - c.x0)) * width}px`;
      band.title = `${options.label ?? "Background"} chirps [${start}, ${stop}); drag to move${options.resizable === false ? "" : ", drag edges to resize"}`;
      if (band.querySelector(".start"))
        band.querySelector(".start").hidden = start < c.x0;
      if (band.querySelector(".stop"))
        band.querySelector(".stop").hidden = stop > c.x1;
      for (const button of band.querySelectorAll("button"))
        button.disabled = !options.enabled();
    }
  }
  set(values, width, height, config) {
    const gl = this.gl;
    if (
      width > gl.getParameter(gl.MAX_TEXTURE_SIZE) ||
      height > gl.getParameter(gl.MAX_TEXTURE_SIZE)
    )
      throw Error(
        "This view exceeds the GPU texture limit. Shorten the chirp interval.",
      );
    this.baseConfig = {...config};
    this.config = {...config,...this.viewport};
    this.updateResetButton();
    this.values = values;
    this.width = width;
    this.height = height;
    gl.bindTexture(gl.TEXTURE_2D, this.texture);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.R32F,
      width,
      height,
      0,
      gl.RED,
      gl.FLOAT,
      values,
    );
    this.draw();
  }
  overlay(config) {
    if (this.config) {
      Object.assign(this.baseConfig, config);
      Object.assign(this.config, config);
      this.draw();
    }
  }
  draw() {
    if (!this.config) return;
    this.drawWindow();
    const { clientWidth: W, clientHeight: H } = this.root;
    if (!W || !H) return;
    const gl = this.gl,
      r = this.pixels.getBoundingClientRect();
    this.pixels.width = Math.round(
      r.width * (exportDensity || Math.min(devicePixelRatio || 1, 2)),
    );
    this.pixels.height = Math.round(
      r.height * (exportDensity || Math.min(devicePixelRatio || 1, 2)),
    );
    gl.viewport(0, 0, this.pixels.width, this.pixels.height);
    gl.useProgram(this.program);
    gl.uniform1f(gl.getUniformLocation(this.program, "lo"), this.config.lo);
    gl.uniform1f(gl.getUniformLocation(this.program, "hi"), this.config.hi);
    gl.uniform1i(
      gl.getUniformLocation(this.program, "kind"),
      this.config.kind || 0,
    );
    const c0=this.baseConfig,c1=this.config;
    gl.uniform4f(gl.getUniformLocation(this.program,'sourceRect'),
      (c1.x0-c0.x0)/(c0.x1-c0.x0),(c1.y0-c0.y0)/(c0.y1-c0.y0),
      (c1.x1-c0.x0)/(c0.x1-c0.x0),(c1.y1-c0.y0)/(c0.y1-c0.y0));
    gl.drawArrays(gl.TRIANGLES, 0, 6);
    const ctx = sizeCanvas(this.axes, W, H);
    ctx.clearRect(0, 0, W, H);
    const left = 74,
      top = 37,
      width = W - 176,
      height = H - 100,
      c = this.config;
    ctx.font = "14px Arial, sans-serif";
    ctx.fillStyle = fg;
    ctx.strokeStyle = grid;
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (let j = 0; j <= 4; j++) {
      let v = c.y0 + ((c.y1 - c.y0) * j) / 4,
        y = top + height * (1 - j / 4);
      ctx.fillText(c.yFormat?.(v) ?? v.toFixed(2), left - 8, y);
      ctx.beginPath();
      ctx.moveTo(left, y);
      ctx.lineTo(left + width, y);
      ctx.stroke();
    }
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    const ticks = Math.max(2, Math.min(6, Math.floor(width / 155)));
    for (let j = 0; j <= ticks; j++) {
      const v = c.x0 + ((c.x1 - c.x0) * j) / ticks,
        x = left + (width * j) / ticks;
      ctx.fillText(
        c.xFormat?.(v) ?? Math.round(v).toString(),
        x,
        top + height + 8,
      );
      if (c.topFormat) ctx.fillText(c.topFormat(v), x, 9);
    }
    ctx.fillText(c.xLabel ?? "", left + width / 2, H - 20);
    ctx.save();
    ctx.translate(8, top + height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(c.yLabel ?? "", 0, 0);
    ctx.restore();
    if (c.topLabel) {
      ctx.textAlign = "right";
      ctx.fillText(c.topLabel, left + width, 24);
    }
    const xp = (v) => left + ((v - c.x0) / (c.x1 - c.x0)) * width;
    for (const band of c.bands ?? []) {
      let a = Math.max(c.x0, band.start),
        b = Math.min(c.x1, band.stop);
      if (b <= a) continue;
      ctx.fillStyle = band.color;
      ctx.fillRect(
        xp(a),
        top,
        Math.max(band.strong ? 2 : 0, xp(b) - xp(a)),
        height,
      );
      ctx.strokeStyle = band.line ?? band.color;
      ctx.lineWidth = band.strong ? 2 : 1;
      ctx.beginPath();
      ctx.moveTo(xp(a), top);
      ctx.lineTo(xp(a), top + height);
      ctx.moveTo(xp(b), top);
      ctx.lineTo(xp(b), top + height);
      ctx.stroke();
      if (band.label) {
        const labelWidth = ctx.measureText(band.label).width + 12;
        const center = Math.max(
          left + labelWidth / 2,
          Math.min(left + width - labelWidth / 2, (xp(a) + xp(b)) / 2),
        );
        ctx.fillStyle = "#0c1923e8";
        ctx.fillRect(
          center - labelWidth / 2,
          top + height - 25,
          labelWidth,
          20,
        );
        ctx.fillStyle = band.line;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(band.label, center, top + height - 15);
      }
    }
    if (c.frames) {
      ctx.strokeStyle = "#b6c9d344";
      ctx.lineWidth = 0.6;
      for (
        let k = Math.ceil(c.x0 / c.frames) * c.frames;
        k < c.x1;
        k += c.frames
      ) {
        ctx.beginPath();
        ctx.moveTo(xp(k), top);
        ctx.lineTo(xp(k), top + height);
        ctx.stroke();
      }
    }
    if (c.marker) {
      const source=this.baseConfig;
      const markerX=source.x0+(c.marker.column+0.5)/this.width*(source.x1-source.x0);
      const markerY=source.y0+(c.marker.row+0.5)/this.height*(source.y1-source.y0);
      const mx=xp(markerX);
      const my=top+(1-(markerY-c.y0)/(c.y1-c.y0))*height;
      ctx.save();
      ctx.beginPath();
      ctx.rect(left, top, width, height);
      ctx.clip();
      for (const [color, lineWidth] of [
        ["#07131c", 5],
        ["#ffffff", 2],
      ]) {
        ctx.strokeStyle = color;
        ctx.lineWidth = lineWidth;
        ctx.beginPath();
        ctx.moveTo(mx - 8, my);
        ctx.lineTo(mx + 8, my);
        ctx.moveTo(mx, my - 8);
        ctx.lineTo(mx, my + 8);
        ctx.stroke();
      }
      ctx.restore();
    }
    const bx = left + width + 15,
      bw = 12;
    for (let i = 0; i < height; i++) {
      const rgb = (c.kind ? magma : diverging)(1 - i / height);
      ctx.fillStyle = `rgb(${rgb.join(",")})`;
      ctx.fillRect(bx, top + i, bw, 1.2);
    }
    ctx.fillStyle = fg;
    ctx.textAlign = "left";
    for (let i = 0; i <= 4; i++)
      ctx.fillText(
        (c.hi - ((c.hi - c.lo) * i) / 4).toFixed(
          c.colorDecimals ?? (c.kind ? 1 : 0),
        ),
        bx + bw + 5,
        top + (height * i) / 4,
      );
    ctx.save();
    ctx.translate(W - 13, top + height / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(c.colorLabel ?? "", 0, 0);
    ctx.restore();
  }
}
export class LinePlot {
  constructor(id, onSelect) {
    this.root = document.getElementById(id);
    this.canvas = document.createElement("canvas");
    this.root.append(this.canvas);
    addExport(this,false);
    this.root.style.display = "none";
    this.canvas.onclick = (e) => {
      if (!this.config || !onSelect) return;
      const r = this.canvas.getBoundingClientRect(),
        w = r.width - 100;
      onSelect(
        this.config.x0 +
          ((e.clientX - r.left - 65) / w) * (this.config.x1 - this.config.x0),
      );
    };
    new ResizeObserver(() => this.draw()).observe(this.root);
  }
  set(series, config) {
    this.series = series;
    this.config = config;
    this.root.style.display = "block";
    this.draw();
  }
  draw() {
    if (!this.config) return;
    const W = this.root.clientWidth,
      H = this.root.clientHeight;
    if (!W || !H) return;
    const c = this.config,
      ctx = sizeCanvas(this.canvas, W, H),
      l = 65,
      t = 28,
      w = W - 100,
      h = H - 83;
    ctx.clearRect(0, 0, W, H);
    let values = this.series.flatMap((s) =>
        Array.from(s.y).filter(Number.isFinite),
      ),
      y0 = c.y0 ?? Math.min(...values),
      y1 = c.y1 ?? Math.max(...values);
    if (y0 === y1) {
      y0--;
      y1++;
    }
    const pad = c.y0 === undefined ? (y1 - y0) * 0.08 : 0;
    y0 -= pad;
    y1 += pad;
    const xp = (x) => l + ((x - c.x0) / (c.x1 - c.x0)) * w,
      yp = (y) => t + h - ((y - y0) / (y1 - y0)) * h;
    ctx.font = "14px Arial, sans-serif";
    ctx.fillStyle = fg;
    ctx.strokeStyle = grid;
    ctx.textAlign = "right";
    for (let j = 0; j <= 4; j++) {
      const y = y0 + ((y1 - y0) * j) / 4;
      ctx.fillText(y.toFixed(c.decimals ?? 0), l - 8, yp(y) + 3);
      ctx.beginPath();
      ctx.moveTo(l, yp(y));
      ctx.lineTo(l + w, yp(y));
      ctx.stroke();
    }
    ctx.textAlign = "center";
    for (let i = 0; i <= 4; i++) {
      const x = c.x0 + ((c.x1 - c.x0) * i) / 4;
      ctx.fillText(c.xFormat?.(x) ?? x.toFixed(2), xp(x), t + h + 20);
    }
    ctx.fillText(c.xLabel ?? "", l + w / 2, H - 8);
    ctx.save();
    ctx.translate(13, t + h / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(c.yLabel ?? "", 0, 0);
    ctx.restore();
    let leg = l;
    for (const s of this.series) {
      ctx.fillStyle = s.color;
      ctx.textAlign = "left";
      ctx.fillText(s.name, leg, 14);
      leg += ctx.measureText(s.name).width + 22;
      ctx.strokeStyle = s.color;
      ctx.lineWidth = s.width ?? 1.6;
      ctx.save();
      ctx.beginPath();
      ctx.rect(l, t, w, h);
      ctx.clip();
      ctx.beginPath();
      let pen = false;
      for (let i = 0; i < s.y.length; i++) {
        if (!Number.isFinite(s.y[i])) {
          pen = false;
          continue;
        }
        const x = s.x[i],
          y = s.y[i];
        if (!pen) {
          ctx.moveTo(xp(x), yp(y));
          pen = true;
        } else ctx.lineTo(xp(x), yp(y));
      }
      ctx.stroke();
      ctx.restore();
    }
  }
}
