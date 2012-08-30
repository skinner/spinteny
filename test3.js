goog.require("goog.asserts");
goog.require("goog.dom");
goog.require("goog.style");
goog.require("goog.debug.Logger");
goog.require("goog.vec.Mat4");
goog.require("goog.vec.Vec3");

var vertShader = 
    [
	"uniform    mat4    orgTransforms[8];",
	"uniform    mat4    cameraInverse;",
	"uniform    mat4    cameraProjection;",

	"attribute  vec3    vertex;",
	"attribute  vec3    normal;",
	"attribute  float   org;",

	//"varying    float   light;",

	"void main(void) {",
	//"  light = dot( normalize( mat3( transform[0].xyz, transform[1].xyz, transform[2].xyz ) * normals ), vec3( 0.0, 0.0, -1.0 ));",
	"  mat4 transform = orgTransforms[int(org)];",
	"  gl_Position = cameraProjection * cameraInverse * transform * vec4( vertex, 1.0 );",
	"}"
    ].join( "\n" );

var fragShader =
    [
	"#ifdef GL_ES",
	"precision highp float;",
	"#endif",		
	
	//"varying float light;",

	"void main( void ) {",
	//"gl_FragColor = texture2D( texture, uv ) * light;",
	"  gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);",
	"}"
    ].join( "\n" );

function Spinteny(container) {
    var container = goog.dom.getElement(container);
    var containerSize = goog.style.getSize(container);

    this.context = new GLOW.Context({
	width: containerSize.width,
	height: containerSize.height
    });

    this.context.setupClear( { red: 1, green: 1, blue: 1 } );
    container.appendChild(this.context.domElement);

    this.genomeHeight = 40;
    this.genomeRadius = 300;
    this.chromSpacing = Math.PI / 18;
    
    this.maxOrgs = 8;
    this.orgTransformFlat = new Float32Array(16 * this.maxOrgs);
    this.orgTransforms = [];
    for (var i = 0; i < this.maxOrgs; i++) {
	this.orgTransforms[i] =
	    this.orgTransformFlat.subarray(i * 16, (i * 16) + 16);
    }
    
    goog.vec.Mat4.makeIdentity(this.orgTransforms[0]);
    goog.vec.Mat4.makeIdentity(this.orgTransforms[1]);
    goog.vec.Mat4.translate(this.orgTransforms[0], 0, 70, 0);
    goog.vec.Mat4.translate(this.orgTransforms[1], 0, -30, 0);
    goog.vec.Mat4.rotateY(this.orgTransforms[1], Math.PI/18);

    var thisObj = this;
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
		    new goog.vec.Vec3.createFromValues(0.0,
						       -thisObj.genomeHeight, 
						       0.0),
		    new goog.vec.Vec3.createFromValues(0.0, 0.0, thisObj.genomeRadius),
		    thisObj.chromSpacing
		);
	    }
	);
    
    var LCBs = [0, 1, 2, 3, 4].map(
	    function(chrId) {
		return [
		    [
			[0, chrId, 0, 4000],
			[1, chrId, 0, 4000]
			//[2, chrId, 0, 4000],
			//[3, chrId, 0, 4000],
		    ],
		    [
			[0, chrId, 6000, 10000],
			[1, chrId, 6000, 10000]
			//[2, chrId, 6000, 10000],
			//[3, chrId, 6000, 10000],
		    ]
		];
	    }
	);
    var joined = LCBs.reduce(function(a, b) {
	return a.concat(b);
    }, []);

    var synVerts = this.LCBsToVertices(joined);

    var anchorShaderInfo = {
	vertexShader: vertShader,
	fragmentShader: fragShader,

	data: {
	    // create uniform data

	    orgTransforms: {
		value: this.orgTransformFlat
	    },
	    cameraInverse: GLOW.defaultCamera.inverse,
	    cameraProjection: GLOW.defaultCamera.projection,

	    // create attribute data

	    vertex: synVerts.anchors.vertex,
	    normal: synVerts.anchors.normal,
	    org: synVerts.anchors.org
	},
	// create element data
	primitives: GL.LINE_STRIP
    };

    var anchors = new GLOW.Shader(anchorShaderInfo);

    // Update the default camera position
    GLOW.defaultCamera.localMatrix.setPosition( 0, 0, 500 );
    GLOW.defaultCamera.update();

    this.context.cache.clear();
    this.context.clear();
    anchors.draw();
}

/**
 * add data for two triangles (6 vertices) to dst starting at offset
 * vertices: array of 4 vec3's, in the order
 *           [top left, top right, bottom left, bottom right]
 * dst: flat array of vertex positions (e.g., to pass to gl.drawArrays
 *      with GL.TRIANGLES)
 */
function trianglesForQuad(vertices, dst, offset) {
    // triangles wind counterclockwise (if you're looking at them from
    // outside the mesh), which is what calcArrayFaceNormals expects
    dst.set(vertices[0], offset + 0 ); // tri 1: top left
    dst.set(vertices[2], offset + 3 ); // tri 1: bottom left
    dst.set(vertices[1], offset + 6 ); // tri 1: top right
    dst.set(vertices[1], offset + 9 ); // tri 2: top right
    dst.set(vertices[2], offset + 12); // tri 2: bottom left
    dst.set(vertices[3], offset + 15); // tri 2: bottom right
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
    var result = goog.vec.Vec3.create();
    for (var i = 0; i < src.length; i += 9) {
	// doing some math myself here rather than using a lib
	// function in order to avoid allocating
	v12[0] = src[i + 3] - src[i + 0];
	v12[1] = src[i + 4] - src[i + 1];
	v12[2] = src[i + 5] - src[i + 2];
	v13[0] = src[i + 6] - src[i + 0];
	v13[1] = src[i + 7] - src[i + 1];
	v13[2] = src[i + 8] - src[i + 2];
	goog.vec.Vec3.cross(v12, v13, result);
	goog.vec.Vec3.normalize(result, result);
	dst.set(result, i + 0);
	dst.set(result, i + 3);
	dst.set(result, i + 6);
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
 *
 * the current convention is that each org's drawing area has
 * top at 0 and bottom at 1, which will get transformed by that org's
 * transformation matrix in the vertex shader
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
        otherOrg:  new Float32Array(3 * twistVertCount),
        org:       new Float32Array(twistVertCount)
    };

    var curAnchor = 0;
    var curTwist = 0;
    var topStart, topEnd, botStart, botEnd;
    var orgId = 0, chrId = 1, start = 2, end = 3;
    var sideVec = goog.vec.Vec3.create();

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
	    twists.org[curTwist * 6 + 1] = org;     // bottom left
	    twists.org[curTwist * 6 + 2] = prevOrg; // top right
	    twists.org[curTwist * 6 + 3] = prevOrg; // top right
	    twists.org[curTwist * 6 + 4] = org;     // bottom left
	    twists.org[curTwist * 6 + 5] = org;     // bottom right
	    // and the twists.otherOrg values are the inverse of
	    // the twists.org values
	    twists.otherOrg[curTwist * 6 + 0] = org;
	    twists.otherOrg[curTwist * 6 + 1] = prevOrg;
	    twists.otherOrg[curTwist * 6 + 2] = org;
	    twists.otherOrg[curTwist * 6 + 3] = org;
	    twists.otherOrg[curTwist * 6 + 4] = prevOrg;
	    twists.otherOrg[curTwist * 6 + 5] = prevOrg;

	    goog.vec.Vec3.subtract(twistVerts[0], twistVerts[1], sideVec);
	    twists.sideVec.set(sideVec, curTwist * 18 + 0); // top left
	    twists.sideVec.set(sideVec, curTwist * 18 + 6); // top right
	    twists.sideVec.set(sideVec, curTwist * 18 + 9); // top right

	    twists.otherVert.set(twistVerts[2], curTwist * 18 + 0);
	    twists.otherVert.set(twistVerts[3], curTwist * 18 + 6);
	    twists.otherVert.set(twistVerts[3], curTwist * 18 + 9);

	    goog.vec.Vec3.subtract(twistVerts[2], twistVerts[3], sideVec);
	    twists.sideVec.set(sideVec, curTwist * 18 + 3);  // bottom left
	    twists.sideVec.set(sideVec, curTwist * 18 + 12); // bottom left
	    twists.sideVec.set(sideVec, curTwist * 18 + 15); // bottom right

	    twists.otherVert.set(twistVerts[0], curTwist * 18 + 3);
	    twists.otherVert.set(twistVerts[0], curTwist * 18 + 12);
	    twists.otherVert.set(twistVerts[1], curTwist * 18 + 15);

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
    this.rotAxis = goog.vec.Vec3.normalize(axis,
					   goog.vec.Vec3.createFloat32());
    // the rotations done by this mapper are left-handed, i.e., 
    // in toSpatial, positive distance is downward, and
    // positive rotations go counterclockwise if you're looking
    // at the circle from above
    goog.vec.Vec3.negate(this.rotAxis, this.rotAxis);
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
    // scale the angle down to leave room for padding
    angle *= ( (2 * Math.PI)
	       / ( (2 * Math.PI) 
		   + ( this.padding * ( this.chroms.length ) ) ) );
    var rotM = goog.vec.Mat4.makeRotate(goog.vec.Mat4.createFloat32(),
					angle,
					this.rotAxis[0],
					this.rotAxis[1],
					this.rotAxis[2]);
    var result = goog.vec.Vec3.createFloat32();
    goog.vec.Mat4.multVec3NoTranslate(rotM, this.origin, result);

    if (0 != distance) {
	var distVec = goog.vec.Vec3.createFloat32();
	goog.vec.Vec3.scale(this.axis, distance, distVec);
	goog.vec.Vec3.add(result, distVec, result);
    }

    return result;
};

//TODO: CylMapper.prototype.fromSpatial = function(spatialPos){}
