const fortran = require("../lib")
    , abs = require("abs")
    ;

const input = process.argv[2];
const param = "-o";
const output = "BUMPER";
if (!input) {
    return console.log("Usage: node others.js file.f");
}

fortran(abs(input), param, output, (err, data) => {
    console.log(err || data);
});
