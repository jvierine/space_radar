import assert from 'node:assert/strict';
globalThis.location={href:'https://juha.no/fmcw/lab/'};
const {initialState,defaultState,decodeState,encodeState}=await import('../web/lab/gui-state.mjs');
assert.deepEqual(initialState,defaultState);assert.equal(initialState.controls.chirp,'3821');assert.equal(initialState.controls.vMax,'600');
assert.deepEqual(initialState.background,{start:3633,stop:3692});assert.deepEqual(initialState.analysis,{start:3704,stop:3880});
assert.deepEqual(initialState.view,{start:3369,stop:3957});assert.equal(initialState.fit,true);
const shared={version:1,controls:{pulses:'16'}};assert.deepEqual(decodeState(encodeState(location.href,shared)),shared);
console.log('PASS: requested default GUI state; explicit shared state overrides it');
