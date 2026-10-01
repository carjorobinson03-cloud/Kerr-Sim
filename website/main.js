//Shader Port to WebGL, I believe an src bridge should work to the html file

import { FRAG_SRC } from './shader.js';
import GUI from 'https://cdn.jsdelivr.net/npm/lil-gui@0.19/+esm';

const VERT_SRC = `#version 300 es
void main() {
  vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
  gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}`;

const canvas = document.getElementById('view');

//lad fix test
const RENDER_WIDTH = 1260;
const RENDER_HEIGHT = Math.round(RENDER_WIDTH * 9 /16);
function resize() {
    //const dpr = Math.min(window.devicePixelRatio, 2);
    //canvas.width  = canvas.clientWidth  * dpr;
    //canvas.height = canvas.clientHeight * dpr;
    canvas.width = RENDER_WIDTH;
    canvas.height = RENDER_HEIGHT;
}
new ResizeObserver(resize).observe(canvas);
resize();
const gl = canvas.getContext('webgl2');

const starfield = new Float32Array(await fetch('starfield2.bin').then(r => r.arrayBuffer()));
const bbtorgb = new Float32Array(await fetch('bbtorgb.bin').then(r => r.arrayBuffer()));
const disctemp = new Float32Array(await fetch('temp_2d.bin').then(r => r.arrayBuffer()));


function compile(source, shaderType) {
    const shader = gl.createShader(shaderType);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
    return shader;
}
function uploadTexture(gl, unit, data, internalFmt, w, h) {
    const fmt = (internalFmt === gl.R32F) ? gl.RED : gl.RGB;
    const tex = gl.createTexture();
    gl.activeTexture(unit);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, internalFmt, w, h, 0, fmt, gl.FLOAT, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    return tex
}

function rISCO(a) {
    const Z1 = 1.0 + Math.cbrt(1.0 - a*a) * (Math.cbrt(1.0 + a) + Math.cbrt(1.0 - a));
    const Z2 = Math.sqrt(3.0*a*a + Z1*Z1);
    return 3.0 + Z2 - Math.sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0*Z2));
}
//do not need r_ISCO_prograde because all the python functions are being used to construct binary bake LUTs

const program = gl.createProgram();
gl.attachShader(program, compile(VERT_SRC, gl.VERTEX_SHADER)); //Possible error source from mismatching arg calls.
gl.attachShader(program, compile(FRAG_SRC, gl.FRAGMENT_SHADER));
gl.linkProgram(program);
if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
gl.useProgram(program);
const U = name => gl.getUniformLocation(program, name);

uploadTexture(gl, gl.TEXTURE0, starfield, gl.RGB32F, 4096, 2048);
gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
uploadTexture(gl, gl.TEXTURE1, bbtorgb,   gl.RGB32F, 1024, 1);
uploadTexture(gl, gl.TEXTURE2, disctemp,  gl.R32F, 1024, 64);

gl.uniform1i(U('starfield'), 0); 
gl.uniform1i(U('bbColor'), 1);
gl.uniform1i(U('discTemp'), 2);
gl.uniform1f(U('M'), 1.0);
gl.uniform1f(U('tanHalfFov'), Math.tan(0.5 * 40.0 * Math.PI / 180.0)); // may as well bring in here, never changes
gl.uniform1f(U('T_LUT_MIN'), 1000.0);
gl.uniform1f(U('T_LUT_MAX'), 20000.0);
gl.uniform1f(U('exposure'), 1.2);
gl.uniform2f(U('aRange'), 0.0, 0.998); //Clipped like this to avoid LUT binary file overloading.

//GUI init

let phi_camera = Math.PI / 6, theta_camera = 85*Math.PI / 180, r_cam = 50.0; //Matching Python side.
const physParams = { a: 0.0, T_peak: 4000 };

canvas.addEventListener('pointermove', e => {
    if (e.buttons !== 1) return;
    phi_camera += e.movementX * 0.005;
    theta_camera = Math.max(0.05, Math.min(Math.PI - 0.05, theta_camera - e.movementY * 0.005));
});

canvas.addEventListener('wheel', e=> {
    r_cam = Math.max(8, Math.min(50, r_cam + e.deltaY * 0.05));
});

const gui = new GUI();
gui.add(physParams,'a', 0, 0.998).name('Spin (a)');
gui.add(physParams,'T_peak', 1000, 20000).name('Peak temp (K)');

function main() {
    const sinT = Math.sin(theta_camera), cosT = Math.cos(theta_camera);
    const sinP = Math.sin(phi_camera), cosP = Math.cos(phi_camera);

    canvas.width = canvas.clientWidth;
    canvas.height = canvas.clientHeight;
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clear(gl.COLOR_BUFFER_BIT);

    gl.uniform1f(U('cam_x'), r_cam*sinT*cosP - physParams.a*sinT*sinP); //dynamic uniforms go in render loop. 
    gl.uniform1f(U('cam_y'), r_cam*sinT*sinP + physParams.a*sinT*cosP);
    gl.uniform1f(U('cam_z'), r_cam*cosT);
    gl.uniform1f(U('a'), physParams.a);
    gl.uniform1f(U('T_peak'), physParams.T_peak);
    gl.uniform1f(U('DISC_IN'), rISCO(physParams.a));
    gl.uniform2f(U('resolution'), canvas.width, canvas.height);

    gl.drawArrays(gl.TRIANGLES, 0, 3);
    requestAnimationFrame(main);
}


requestAnimationFrame(main);


