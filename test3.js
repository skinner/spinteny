goog.require("goog.asserts");
goog.require("goog.dom");
goog.require("goog.debug.Logger");
goog.require("goog.vec.Mat4");
goog.require("goog.vec.Vec3");

function Spinteny(container) {
    var container = goog.dom.getElement(container);
    this.context = new GLOW.Context();
    this.context.setupClear( { red: 1, green: 1, blue: 1 } );
    container.appendChild(this.context.domElement);

    this.anchorHeight = 10;
    this.genomeSpacing = 5;

    this.mappers = 
	[0, 1, 2, 3].map(
	    function(orgId) {
		return new CylMapper(
		    [
			{ start: 0, end: 10000, name: "a" },
			{ start: 0, end: 10000, name: "b" },
			{ start: 0, end: 10000, name: "c" },
			{ start: 0, end: 10000, name: "d" },
			{ start: 0, end: 10000, name: "e" }
		    ],
		    new goog.vec.Vec3.createFromValues(0.0, -1.0, 0.0),
		    new goog.vec.Vec3.createFromValues(0.0, 0.0, -1.0)
		);
	    }
	);

    var LCBs = [0, 1, 2, 3, 4].map(
	    function(chrId) {
		return [
		    [
			[0, chrId, 0, 4000],
			[1, chrId, 0, 4000],
			[2, chrId, 0, 4000],
			[3, chrId, 0, 4000],
		    ],
		    [
			[0, chrId, 6000, 10000],
			[1, chrId, 6000, 10000],
			[2, chrId, 6000, 10000],
			[3, chrId, 6000, 10000],
		    ]
		];
	    }
	);
    var joined = LCBs.reduce(function(a, b) {
	return a.concat(b);
    }, []);

    var synVerts = this.LCBsToVertices(joined);
    
    console.log(synVerts);
}

function copyVec3(src, dst, offset) {
    dst[offset + 0] = src[0];
    dst[offset + 1] = src[1];
    dst[offset + 2] = src[2];
}

/**
 * add data for two triangles (6 vertices) to dst starting at offset
 * vertices: array of 6 vec3's
 * dst: flat array of vertex positions (e.g., to pass to gl.drawArrays
 *      with GL.TRIANGLES)
 */
function trianglesForQuad(vertices, dst, offset) {
    copyVec3(vertices[0], dst, offset + 0 ); // tri 1: top left
    copyVec3(vertices[1], dst, offset + 3 ); // tri 1: top right
    copyVec3(vertices[2], dst, offset + 6 ); // tri 1: bottom left
    copyVec3(vertices[1], dst, offset + 9 ); // tri 2: top right
    copyVec3(vertices[2], dst, offset + 12); // tri 2: bottom left
    copyVec3(vertices[3], dst, offset + 15); // tri 2: bottom right
}

/**
 * calculate normals for the vertices in src
 * (assumes src and dst are arrays that could be passed to
 *  gl.drawArrays with GL_TRIANGLES)
 * (also assumes that normals should be perpendicular to the triangle face)
 */
function calcArrayFaceNormals(src, dst) {
    goog.asserts.assert(src.length == dst.length);

    var v12 = goog.vec.Vec3.create();
    var v13 = goog.vec.Vec3.create();
    for (var i = 0; i < src.length; i += 9) {
	v12[0] = src[i + 3] - src[i + 0];
	v12[1] = src[i + 4] - src[i + 1];
	v12[2] = src[i + 5] - src[i + 2];
	v13[0] = src[i + 6] - src[i + 0];
	v13[1] = src[i + 7] - src[i + 1];
	v13[2] = src[i + 8] - src[i + 2];
	var result = goog.vec.Vec3.create();
	goog.vec.Vec3.cross(v12, v13, result);
	goog.vec.Vec3.normalize(result, result);
	copyVec3(result, dst, i + 0);
	copyVec3(result, dst, i + 3);
	copyVec3(result, dst, i + 6);
    }
}

function repeat(val, dst, offset, count) {
    for (var i = 0; i < count; i++) dst[offset + i] = val;
}

/**
 * takes a single genome match in the form of an array:
 * [orgID, chrID, start, end]
 * and returns a quadrilateral that covers that region of the
 * spinteny cylinder in the form of an array of four vec3's:
 * [topStart, topEnd, botStart, botEnd]
 */
Spinteny.prototype.matchToQuad = function(match) {
    return [
	this.mappers[match[0]].toSpatial(match[1], match[2], 0),
	this.mappers[match[0]].toSpatial(match[1], match[3], 0),
	this.mappers[match[0]].toSpatial(match[1], match[2], 1),
	this.mappers[match[0]].toSpatial(match[1], match[3], 1)
    ];
}

/**
 * takes an array of locally collinear blocks of the format:
 * [ [orgID, chrID, start, end], [orgID, chrID, start, end], ... ]
 *
 * and returns an object with
 * { 
 *     anchors:
 *         {
 *             vertex: 
 *             normal:
 *             org:
 *         }
 *     twists:
 *         {
 *             vertex:
 *             org:
 *             sideVec:
 *             otherVert:
 *             otherOrg:
 *         }
 * }
 *
 * "anchors" visually represent the span of an LCB on a particular genome.
 * "twists" visually represent the connections between related anchors.
 * "twists" has extra vertex attributes because we need to calculate
 * vertex normals within the vertex shader for the twists
 * (the normals depend on the position of a twist's top and bottom anchors,
 * and those positions are calculated within the vertex shader)
 * the normals for the twist vertices will be the cross product of
 * sideVec and vector from the current vertex to the otherVert, where
 * the otherVert has been transformed according to the transformation for
 * otherOrg.
 */
Spinteny.prototype.LCBsToVertices = function(blocks) {
    var numAnchors = 0;
    var numTwists = 0;
    for (var i = 0; i < blocks.length; i++) {
	// there's an anchor for each block,
	numAnchors += blocks[i].length;
	// and there's a twist between related anchors
	numTwists += blocks[i].length - 1;
    }
    // anchors and twists are planar quadrilaterals
    // using drawArrays means 6 vertices per face
    var anchVertCount = 6 * numAnchors;
    var anchors = {
	vertex: new Float32Array(3 * anchVertCount),
	normal: new Float32Array(3 * anchVertCount),
	// ideally, org would be an unsigned int (or short or byte) array,
	// but GLSL (in webGL 1.0) doesn't allow those types for attributes
	org:    new Float32Array(anchVertCount)
    };

    var twistVertCount = 6 * numTwists;
    var twists = {
	vertex:    new Float32Array(3 * twistVertCount),
        sideVec:   new Float32Array(3 * twistVertCount),
        otherVert: new Float32Array(3 * twistVertCount),
        otherOrg:  new Float32Array(3 * twistVertCount)
        org:       new Float32Array(twistVertCount),
    };

    var curAnchor = 0;
    var curTwist = 0;
    var topStart, topEnd, botStart, botEnd;
    var orgId = 0, chrId = 1, start = 2, end = 3; 

    for (var blockIdx = 0; blockIdx < blocks.length; blockIdx++) {
	var block = blocks[blockIdx];
	var match = 0;

	// take the genome-space match and convert it into
	// four 3d-space vertex positions for the quadrilateral
	// that will visually represent the match
	var anchorVerts = this.matchToQuad(block[match]);

	// add two triangles for this quad to anchors.vertex
	// (18 is for 3 values (x, y, z) for 6 vertices)
	trianglesForQuad(anchorVerts, anchors.vertex, curAnchor * 18);
	// set anchors.org to this organism ID for all 6 vertices
	// in the quad
	var org = block[match][orgId];
	repeat(org, anchors.org, curAnchor * 6, 6);
	curAnchor++;
	
	for (match = 1; match < block.length; match++) {
	    var oldAnchorVerts = anchorVerts;
	    anchorVerts = this.matchToQuad(block[match]);

	    // the top two vertices for twistVerts are the bottom two
	    // from the previous anchor, and the bottom two vertices for
	    // twistVerts are the top two vertices from the next anchor
	    //
	    //       0            1
	    //       |  oldAnchor |
	    //       2____________3
	    //       /            /
	    //      /    twist   /
	    //     /___________ /
	    //    0            1
            //    |   anchor   |
	    //    2            3
	    twistVerts = [oldAnchorVerts[2], oldAnchorVerts[3],
			  anchorVerts[0], anchorVerts[1]];
	    trianglesForQuad(twistVerts, twists.vertex, curTwist * 18);

	    var prevOrg = org;
	    org = block[match][orgId];
	    // the top vertices of the twist are for the previous org,
	    // and the bottom vertices of the twist are for the current org
	    twists.org[curTwist * 6 + 0] = prevOrg; // top left
	    twists.org[curTwist * 6 + 1] = prevOrg; // top right
	    twists.org[curTwist * 6 + 2] = org;     // bottom left
	    twists.org[curTwist * 6 + 3] = prevOrg; // top right
	    twists.org[curTwist * 6 + 4] = org;     // bottom left
	    twists.org[curTwist * 6 + 5] = org;     // bottom right
	    // and the twists.otherOrg values are the inverse of
	    // the twists.org values
	    twists.otherOrg[curTwist * 6 + 0] = org;
	    twists.otherOrg[curTwist * 6 + 1] = org;
	    twists.otherOrg[curTwist * 6 + 2] = prevOrg;
	    twists.otherOrg[curTwist * 6 + 3] = org;
	    twists.otherOrg[curTwist * 6 + 4] = prevOrg;
	    twists.otherOrg[curTwist * 6 + 5] = prevOrg;

	    var sideVec = goog.vec.Vec3.create();
	    goog.vec.Vec3.subtract(twistVerts[0], twistVerts[1], sideVec);
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 0); // TL
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 3); // TR
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 9); // TR

	    copyVec3(twistVerts[2], twists.otherVert, curTwist * 18 + 0);
	    copyVec3(twistVerts[3], twists.otherVert, curTwist * 18 + 3);
	    copyVec3(twistVerts[3], twists.otherVert, curTwist * 18 + 9);

	    goog.vec.Vec3.subtract(twistVerts[2], twistVerts[3], sideVec);
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 6);  // BL
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 12); // BL
	    copyVec3(sideVec, twists.sideVec, curTwist * 18 + 15); // BR
	    
	    copyVec3(twistVerts[0], twists.otherVert, curTwist * 18 + 6);
	    copyVec3(twistVerts[0], twists.otherVert, curTwist * 18 + 12);
	    copyVec3(twistVerts[1], twists.otherVert, curTwist * 18 + 15);
	    
	    curTwist++;

	    trianglesForQuad(anchorVerts, anchors.vertex, curAnchor * 18);
	    repeat(org, anchors.org, curAnchor * 6, 6);
	    curAnchor++;
	}
    }

    calcArrayFaceNormals(anchors.vertex, anchors.normal);

    return { anchors: anchors, twists: twists };
}

/**
 * maps between genomic and spatial coordinates on a cylinder
 * @param chroms array of {start, end, name (optional)} objects
 * @param axis {Vec3} axis vector at the center of the cylinder
 * @param origin {Vec3} vector from the axis to the zero point
 * @param padding spacing between chroms (radians) (optional)
 */
function CylMapper(chroms, axis, origin, padding) {
    this.axis = axis;
    this.origin = origin;
    // default padding of 5 degrees
    this.padding = (padding === undefined) ? Math.PI/36 : padding;
    this.chroms = chroms;
    //this.partialSums[i] = total length of chroms[0..(i-1)]
    //this.partialSums[0] = 0
    this.partialSums = [];
    this.byName = {};

    var partialSum = 0;
    for (var i = 0; i < chroms.length; i++) {
	var name = ("name" in chroms[i]) ? chroms[i].name : i;
	this.byName["" + name] = i;
	this.partialSums.push(partialSum);
	partialSum += chroms[i].end - chroms[i].start;
    }
    this.totalLength = partialSum;
}

/**
 * takes a genomic location and converts it to a spatial location
 * @param index index of the chrom in this mapper's chrom list
 * @param base position to convert
 * @param distance (optional) distance along the cylinder's axis
 * @return vec3 with x,y,z of result
 */
CylMapper.prototype.toSpatial = function(index, base, distance) {
    distance = (distance === undefined) ? 0 : distance;
    //var index = this.byName["" + chrom];
    var angle = 
	( ( ( this.partialSums[index] + ( base - this.chroms[index].start ) )
	    / this.totalLength )
	  * ( 2 * Math.PI ) )
	+ (this.padding * index);
    var rotM = goog.vec.Mat4.makeRotate(goog.vec.Mat4.create(), angle,
					this.axis[0],
					this.axis[1],
					this.axis[2]);
    var result = goog.vec.Vec3.create();
    goog.vec.Mat4.multVec3(rotM, this.origin, result);

    if (0 != distance) {
	var distVec = goog.vec.Vec3.create();
	goog.vec.Vec3.scale(this.axis, distance, distVec);
	goog.vec.Vec3.add(result, distVec, result);
    }

    return result;
};

//TODO: CylMapper.prototype.fromSpatial = function(spatialPos){}
