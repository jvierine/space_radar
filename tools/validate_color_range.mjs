import assert from 'node:assert/strict';
import {visibleColorRange} from '../web/lab/color-range.mjs';
const bounds={x0:0,x1:4,y0:0,y1:2};
const data=[0,20,30,100,-Infinity,22,28,NaN];
assert.deepEqual(visibleColorRange(data,4,2,bounds,bounds),{lo:0,hi:100});
assert.deepEqual(visibleColorRange(data,4,2,bounds,{x0:1,x1:3,y0:0,y1:2}),{lo:20,hi:30});
assert.deepEqual(visibleColorRange([5,6,7,8],2,2,{x0:0,x1:2,y0:0,y1:2},{x0:0,x1:2,y0:0,y1:2}),{lo:5,hi:8});
assert.deepEqual(visibleColorRange([0,0],2,1,{x0:0,x1:2,y0:0,y1:1},{x0:0,x1:2,y0:0,y1:1}),{lo:0,hi:1});
console.log('PASS: actual extrema, zoomed extrema, positive minima, missing cells and flat maps');
