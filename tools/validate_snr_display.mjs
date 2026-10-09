import assert from 'node:assert/strict';
import {snrDb,estimateRcs} from '../web/lab/rcs.mjs';
// Display floor is separate from physical power estimation.
for(const ratio of [0,0.5,1,1.5,2]) assert.equal(snrDb(ratio),0);
assert.equal(snrDb(11),10);
assert.equal(snrDb(101),20);
const parameters={temperature:9000,range:1,frequency:77e9,T_coh:144e-6,txPowerDbm:12,txGainDbi:6,rxGainDbi:6,lossDb:0};
assert.equal(estimateRcs(1,parameters),0);
assert.equal(estimateRcs(1.5,parameters)/estimateRcs(2,parameters),0.5);
console.log('PASS: 0 dB display floor; positive SNR and unclamped RCS powers preserved');
