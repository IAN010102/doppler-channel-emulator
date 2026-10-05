'use strict';
// Regenerates tests/golden/default_snapshot.json. Run it ONLY when a change of the numbers is intended, and say so in the commit message.
const fs = require('fs'), path = require('path');
const { compute } = require('./golden_cases.js');
const file = path.join(__dirname, 'golden', 'default_snapshot.json');
fs.writeFileSync(file, JSON.stringify(compute(), null, 1) + '\n');
console.log('written', file);
