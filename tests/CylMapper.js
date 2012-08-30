goog.require("goog.vec.Vec3");

//TODO actually integrate with a testing framework
function vecApproxEq(a, b, tolerance) {
    b = goog.vec.Vec3.createFloat32FromArray(b);
    var dist = goog.vec.Vec3.distance(a, b);
    if (dist > tolerance) console.log([a, "did not match", b]);
    return dist <= tolerance;
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
    vecApproxEq(testMapper.toSpatial(0, 0,     0), [ 0,    0, 100  ], 1e-12),
    vecApproxEq(testMapper.toSpatial(0, 10000, 0), [ 100,  0, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(1, 0,     0), [ 100,  0, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(1, 10000, 0), [ 0,    0, -100 ], 1e-12),
    vecApproxEq(testMapper.toSpatial(2, 0,     0), [ 0,    0, -100 ], 1e-12),
    vecApproxEq(testMapper.toSpatial(2, 10000, 0), [ -100, 0, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(3, 0,     0), [ -100, 0, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(3, 10000, 0), [ 0,    0, 100  ], 1e-12),

    vecApproxEq(testMapper.toSpatial(0, 0,     1), [ 0,    -100, 100  ], 1e-12),
    vecApproxEq(testMapper.toSpatial(0, 10000, 1), [ 100,  -100, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(1, 0,     1), [ 100,  -100, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(1, 10000, 1), [ 0,    -100, -100 ], 1e-12),
    vecApproxEq(testMapper.toSpatial(2, 0,     1), [ 0,    -100, -100 ], 1e-12),
    vecApproxEq(testMapper.toSpatial(2, 10000, 1), [ -100, -100, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(3, 0,     1), [ -100, -100, 0    ], 1e-12),
    vecApproxEq(testMapper.toSpatial(3, 10000, 1), [ 0,    -100, 100  ], 1e-12)
]);
