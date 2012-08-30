//TODO actually integrate with a testing framework
function arrayApproxEq(arr1, arr2, tolerance) {
    tolerance = (tolerance === undefined) ? 1e-12 : tolerance;
    if (arr1.length != arr2.length) {
	console.log("arrays have different length");
	return false
    }
    for (var i = 0; i < arr1.length; i++) {
	if (Math.abs(arr1[i] - arr2[i]) >= tolerance) {
	    console.log("not equal at " + i);
	    return false;
	}
    }
    return true;
}

var testMapper = new CylMapper(
    [
	{ start: 0, end: 10000, name: "a" },
	{ start: 0, end: 10000, name: "b" },
	{ start: 0, end: 10000, name: "c" },
	{ start: 0, end: 10000, name: "d" }
    ],
    new goog.vec.Vec3.createFromValues(0.0, -100.0, 0.0),
    new goog.vec.Vec3.createFromValues(0.0, 0.0, 100.0),
    0.0
);

console.log([
    arrayApproxEq(testMapper.toSpatial(0, 0,     0), [ 0,    0, 100  ]),
    arrayApproxEq(testMapper.toSpatial(0, 10000, 0), [ 100,  0, 0    ]),
    arrayApproxEq(testMapper.toSpatial(1, 0,     0), [ 100,  0, 0    ]),
    arrayApproxEq(testMapper.toSpatial(1, 10000, 0), [ 0,    0, -100 ]),
    arrayApproxEq(testMapper.toSpatial(2, 0,     0), [ 0,    0, -100 ]),
    arrayApproxEq(testMapper.toSpatial(2, 10000, 0), [ -100, 0, 0    ]),
    arrayApproxEq(testMapper.toSpatial(3, 0,     0), [ -100, 0, 0    ]),
    arrayApproxEq(testMapper.toSpatial(3, 10000, 0), [ 0,    0, 100  ]),

    arrayApproxEq(testMapper.toSpatial(0, 0,     1), [ 0,    -100, 100  ]),
    arrayApproxEq(testMapper.toSpatial(0, 10000, 1), [ 100,  -100, 0    ]),
    arrayApproxEq(testMapper.toSpatial(1, 0,     1), [ 100,  -100, 0    ]),
    arrayApproxEq(testMapper.toSpatial(1, 10000, 1), [ 0,    -100, -100 ]),
    arrayApproxEq(testMapper.toSpatial(2, 0,     1), [ 0,    -100, -100 ]),
    arrayApproxEq(testMapper.toSpatial(2, 10000, 1), [ -100, -100, 0    ]),
    arrayApproxEq(testMapper.toSpatial(3, 0,     1), [ -100, -100, 0    ]),
    arrayApproxEq(testMapper.toSpatial(3, 10000, 1), [ 0,    -100, 100  ])
]);
