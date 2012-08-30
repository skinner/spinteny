var tri = new Float32Array([0, 0, 0,
			    1, 0, 0,
			    0, 1, 0]);
var result = new Float32Array(9);
calcArrayFaceNormals(tri, result);
console.log(arrayApproxEq(result, new Float32Array([0, 0, 1,
						    0, 0, 1,
						    0, 0, 1])));
tri = new Float32Array([0, 0, 0,
			0, 1, 0,
			1, 0, 0]);
calcArrayFaceNormals(tri, result);
console.log(arrayApproxEq(result, new Float32Array([0, 0, -1,
						    0, 0, -1,
						    0, 0, -1])));

