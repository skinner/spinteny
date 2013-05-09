/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

goog.require("goog.asserts");
goog.require("goog.dom");
goog.require("goog.debug.Logger");
goog.require("goog.events");
goog.require("goog.style");
goog.require("goog.ui.Component");
goog.require("goog.ui.Slider");
goog.require("goog.vec.Mat4");
goog.require("goog.vec.Vec3");

var anchorVertShader = 
    [
        "uniform  mat4  orgTransforms[8];",
        "uniform  mat4  viewMatrix;",
        "uniform  mat4  projMatrix;",

        "attribute  vec3   vertex;",
        "attribute  vec3   normal;",
        "attribute  float  org;",

        "varying  vec3   vNormalEye;",

        "void main(void) {",
        "  mat4 mvMatrix = viewMatrix * orgTransforms[int(org)];",
        //"  vNormalEye = transpose(inverse(mat3(mvMatrix))) * normal;",
        "  vNormalEye = mat3(mvMatrix) * normal;",
        "  gl_Position = projMatrix * mvMatrix * vec4( vertex, 1.0 );",
        "}"
    ].join( "\n" );

var twistVertShader = 
    [
        "uniform    mat4    orgTransforms[8];",
        "uniform    mat4    viewMatrix;",
        "uniform    mat4    projMatrix;",

        "attribute  vec3    start;",
        "attribute  vec3    end;",
        "attribute  vec3    otherStart;",
        "attribute  vec3    otherEnd;",
        "attribute  float   prevOrg;",
        "attribute  float   nextOrg;",
        "attribute  float   distance;",
        "attribute  float   normMult;",
        "varying vec3 vNormalEye;",

        "void main(void) {",
        // twist vertices are all along the left and right sides of
        // each twist.  thisStart and thisEnd are the start and end
        // points of the line this vertex is on.  otherStart and
        // otherEnd are the start and end points of the opposide side
        // of the twist.  distance ([0, 1]) is how far along those
        // lines this vertex goes.
        //
        // To calculate the normal, we take the cross product of the
        // line we're on and a vector across the twist to the point at
        // the same "distance" along the other edge of the twist.
        // see dev-notes/twist-shader.png

        "  mat4 startTrans = orgTransforms[int(prevOrg)];",
        "  mat4 endTrans = orgTransforms[int(nextOrg)];",
        "  vec4 thisStart = startTrans * vec4(start, 1.0);",
        "  vec4 thisEnd = endTrans * vec4(end, 1.0);",
        "  vec4 thisPos = thisStart + (distance * (thisEnd - thisStart));",

        "  vec4 otherStart = startTrans * vec4(otherStart, 1.0);",
        "  vec4 otherEnd = endTrans * vec4(otherEnd, 1.0);",
        "  vec4 otherPos = otherStart + (distance * (otherEnd - otherStart));",

        "  vec4 thisVec = thisEnd - thisPos;",
        "  vec4 otherVec = otherPos - thisPos;",
        "  vNormalEye = normMult * cross(thisVec.xyz, otherVec.xyz);",

        "  gl_Position = projMatrix * viewMatrix * vec4(thisPos.xyz, 1.0);",
        "}"
    ].join( "\n" );


var alignFragShader =
    [
        "#ifdef GL_ES",
        "precision highp float;",
        "#endif",               

        "uniform  vec3   ambColor;",
        "uniform  vec3   dirColor;",
        "uniform  vec3   dirVector;",
        "uniform  float  ambWeight;",
        "uniform  float  dirWeight;",
        "uniform  vec2   fogRange;",

        "varying  vec3   vNormalEye;",

        "void main() {",
        "    float z = gl_FragCoord.z / gl_FragCoord.w;",
        "    float fogFactor = (fogRange.y - z) / (fogRange.y - fogRange.x);",
        "    fogFactor = clamp(fogFactor, 0.0, 1.0);",

        "    vec3 normal = normalize(vNormalEye);",
        // flip normal for rear-facing fragments, because we want our
        // materials to be double-sided
        "    normal = normal * ( -1.0 + 2.0 * float( gl_FrontFacing ) );",

        "    float dirDiffuse = max(dot(normal, dirVector), 0.0);",
        "    vec3 dirColor = dirWeight * dirDiffuse * dirColor;",
        "    vec3 ambColor = ambWeight * ambColor;",
        "    vec3 color = dirColor + ambColor;",
        "    gl_FragColor = vec4(mix(vec3(1.0), color, fogFactor), 1.0);",
        "}"
    ].join( "\n" );

function Spinteny(container, orgChroms, LCBs) {
    this.container = goog.dom.getElement(container);
    this.containerSize = goog.style.getSize(this.container);

    this.context = new GLOW.Context({
        width: this.containerSize.width,
        height: this.containerSize.height
    });
    this.canvas = this.context.domElement;
    this.canvas.style.width = "100%";
    this.canvas.style.height = "100%";

    this.context.setupClear( { red: 1, green: 1, blue: 1, alpha: 1 } );
    this.container.appendChild(this.canvas);

    this.ambColor = goog.vec.Vec3.createFromValues(0.75, 0.75, 0.85);
    this.dirColor = goog.vec.Vec3.createFromValues(0.70, 0.75, 0.85);
    this.dirVector = goog.vec.Vec3.createFromValues(-0.4, 0.1, 1.0);
    goog.vec.Vec3.normalize(this.dirVector, this.dirVector);
    this.dirWeight = new Float32Array([0.5]);
    this.ambWeight = new Float32Array([0.5]);

    this.nearClip = 10;
    this.farClip = 10000;

    this.cameraDistance = 100; // arbitrary, must be larger than nearClip
    this.aspect = this.containerSize.width / this.containerSize.height;

    this.fovY = 20 * (Math.PI / 180); // in radians
    this.fovX = 2 * Math.atan(Math.tan(this.fovY / 2) * this.aspect);

    this.viewMatrix = goog.vec.Mat4.createFloat32();
    this.projMatrix = goog.vec.Mat4.createFloat32();
    this.updateViewProj()

    this.genomeCount = 3;

    // zoomMultiplier is an experiment to try and work around numeric
    // precision issues at high zoom levels.  Doesn't seem to
    // make a difference, but I'm leaving it in for now in case
    // it's useful for later experimentation
    this.zoomMultiplier = 1;

    this.minZoom = 0.8;
    this.maxZoom = 200000;

    // a genomeRadius of (this.cameraDistance * radiusFactor)
    // completely fills the container horizontally.
    // see dev-notes/genome-radius.jpg
    var radiusFactor = ( 1 / ( (1 / Math.sin(this.fovX / 2)) - 1 ) );
    this.genomeRadius = (this.cameraDistance * radiusFactor) / this.zoomMultiplier;

    this.genomeHeight =
        this.heightAt(this.cameraDistance) / ((this.genomeCount * 2) - 1);

    this.chromSpacing = Math.PI / 18;

    this.zoomFactor = this.zoomMultiplier;

    this.rotations = new Float32Array(this.genomeCount);
    this.zooms = new Float32Array(this.genomeCount);

    this.orgTransformFlat = new Float32Array(16 * this.genomeCount);
    this.orgTransforms = [];
    for (var i = 0; i < this.genomeCount; i++) {
        this.orgTransforms[i] =
            this.orgTransformFlat.subarray(i * 16, (i * 16) + 16);

        this.zooms[i] = this.zoomFactor;
        this.updateTransform(i);
    }
    // just for dummy data
    this.rotations[1] = Math.PI/18;
    this.rotations[2] = Math.PI/1.6;
    this.updateTransform(1);
    this.updateTransform(2);    

    // when calculating the far end of the fogRange, we multiply
    // the distance by fogMultiplier; changing fogMultiplier
    // changes how visible the distant parts of the cylinder are.
    this.fogMultiplier = 1.8;
    // fragments fade from their normal color at fogRange[0] to
    // white at fogRange[1]
    this.fogRange = new Float32Array([
        this.cameraDistance,
        this.cameraDistance + (this.genomeRadius
                               * this.fogMultiplier
                               * this.zoomMultiplier)
    ]);

    var thisObj = this;
    this.mappers = 
        orgChroms.map(
            function(chromData) {
                return new CylMapper(
                    chromData,
                    new goog.vec.Vec3.createFromValues(0.0,
                                                       -thisObj.genomeHeight, 
                                                       0.0),
                    new goog.vec.Vec3.createFromValues(0.0, 0.0, thisObj.genomeRadius),
                    thisObj.chromSpacing
                );
            }
        );

    var synVerts = this.LCBsToVertices(LCBs);

    var anchorShaderInfo = {
        vertexShader: anchorVertShader,
        fragmentShader: alignFragShader,

        data: {
            // uniforms

            orgTransforms: { value: this.orgTransformFlat },
            viewMatrix: { value: this.viewMatrix },
            projMatrix: { value: this.projMatrix },
            fogRange: { value: this.fogRange },

            ambColor: { value: this.ambColor },
            ambWeight: { value: this.ambWeight },
            dirColor: { value: this.dirColor },
            dirWeight: { value: this.dirWeight },
            dirVector: { value: this.dirVector },

            // attributes

            vertex: synVerts.anchors.vertex,
            normal: synVerts.anchors.normal,
            org: synVerts.anchors.org,
        },
        primitives: GL.TRIANGLES
    };

    var twistShaderInfo = {
        vertexShader: twistVertShader,
        fragmentShader: alignFragShader,

        data: {
            // uniforms

            orgTransforms: { value: this.orgTransformFlat },
            viewMatrix: { value: this.viewMatrix },
            projMatrix: { value: this.projMatrix },
            fogRange: { value: this.fogRange },

            ambColor: { value: this.ambColor },
            ambWeight: { value: this.ambWeight },
            dirColor: { value: this.dirColor },
            dirWeight: { value: this.dirWeight },
            dirVector: { value: this.dirVector },

            // attributes

            start: synVerts.twists.start,
            end: synVerts.twists.end,
            otherStart: synVerts.twists.otherStart,
            otherEnd: synVerts.twists.otherEnd,
            prevOrg: synVerts.twists.prevOrg,
            nextOrg: synVerts.twists.nextOrg,
            distance: synVerts.twists.distance,
            normMult: synVerts.twists.normMult,
        },
        primitives: GL.TRIANGLE_STRIP
    };

    this.anchors = new GLOW.Shader(anchorShaderInfo);
    this.twists = new GLOW.Shader(twistShaderInfo);

    GL.disable(GL.CULL_FACE);

    this.drag = {};
    this.setDragHandler();
    this.addUI();
}

Spinteny.prototype.addUI = function() {
    this.ui = {};

    this.ui.zoomSlider = new goog.ui.Slider;
    this.ui.zoomSlider.setOrientation(goog.ui.Slider.Orientation.VERTICAL);

    this.ui.zoomSlider.setUnitIncrement(0.03);
    this.ui.zoomSlider.setStep(null);
    this.ui.zoomSlider.setMoveToPointEnabled(true);
    //var zoomExp = 14;
    //this.ui.zoomSlider.setMinimum(Math.pow(this.minZoom, 1/zoomExp));
    this.ui.zoomSlider.setMinimum(1);
    //this.ui.zoomSlider.setMaximum(Math.pow(this.maxZoom, 1/zoomExp));
    this.ui.zoomSlider.setMaximum(Math.log(this.maxZoom / this.minZoom)
                                  + 1);
    //this.ui.zoomSlider.setValue(Math.pow(this.minZoom, 1/zoomExp));
    this.ui.zoomSlider.setValue(1);

    this.ui.zoomSlider.createDom();
    var el = this.ui.zoomSlider.getElement();
    el.style.position = "absolute";
    el.style.top = "20px";
    el.style.left = "0px";
    el.style.width = "40px";
    el.style.bottom = "20px";
    this.ui.zoomSlider.render(this.container);

    var thisObj = this;
    function updateZoom() {
        thisObj.zoomFactor = 
            Math.exp(thisObj.ui.zoomSlider.getValue())
            * (thisObj.minZoom / Math.E);
        thisObj.zoomFactor *= thisObj.zoomMultiplier;
        //thisObj.fogRange[1] = thisObj.fogRange[0] + ((thisObj.genomeRadius * thisObj.fogMultiplier * thisObj.zoomMultiplier) / thisObj.zoomFactor);
        //thisObj.fogRange[1] = thisObj.fogRange[0] + (thisObj.genomeRadius * thisObj.fogMultiplier * thisObj.zoomMultiplier);

        for (var i = 0; i < thisObj.genomeCount; i++) {
            thisObj.zooms[i] = thisObj.zoomFactor;
            thisObj.updateTransform(i);
        }

        thisObj.draw();
    }
    this.ui.zoomSlider.addEventListener(
        goog.ui.Component.EventType.CHANGE, updateZoom
    );
    updateZoom();
};

Spinteny.prototype.updateTransform = function(i) {
    goog.vec.Mat4.makeIdentity(this.orgTransforms[i]);
    
    goog.vec.Mat4.translate(this.orgTransforms[i],
                            0,
                            (this.heightAt(this.cameraDistance) / 2)
                            - (i * 2 * this.genomeHeight),
                            -(this.genomeRadius * this.zoomMultiplier));
                            //-(this.genomeRadius * this.zoomMultiplier * this.zooms[i])); // this version also zooms in Z
                            //-(this.genomeRadius * this.zoomMultiplier * (1 + ((this.zooms[i] - 1) * 0.01)))); // this version also zooms a little bit in Z

    goog.vec.Mat4.scale(this.orgTransforms[i],
                        this.zooms[i],
                        1,
                        this.zoomMultiplier);
                        //this.zoomMultiplier * this.zooms[i]); // this version also zooms in Z
                        //this.zoomMultiplier * (1 + ((this.zooms[i] - 1) * 0.01))); // this version also zooms a little bit in Z

    goog.vec.Mat4.rotateY(this.orgTransforms[i], this.rotations[i]);
};

Spinteny.prototype.draw = function() {
    this.context.cache.clear();
    this.context.clear();
    this.anchors.draw();
    this.twists.draw();
};

Spinteny.prototype.updateViewProj = function() {
    goog.vec.Mat4.makeLookAt(
        this.viewMatrix,
        [0, 0, this.cameraDistance], // eye point
        [0, 0, -1],                  // center point
        [0, 1, 0]                    // world up vector
    );

    goog.vec.Mat4.makePerspective(
        this.projMatrix,
        this.fovY,     // y-axis FOV in radians
        this.aspect,   // aspect ratio (width/height)
        this.nearClip, // distance to near clipping plane
        this.farClip   // distance to far clipping plane
    );
};

Spinteny.prototype.heightAt = function(distance) {
    return 2 * distance * Math.tan(this.fovY / 2);
};

Spinteny.prototype.widthAt = function(distance) {
    return this.aspect * this.heightAt(distance);
};

//x and y on [-1, 1]
Spinteny.prototype.coordAtWorldZ = function(x, y, worldZ) {
    var pvInvM = goog.vec.Mat4.create();
    goog.vec.Mat4.multMat(this.projMatrix, this.viewMatrix, pvInvM);
    goog.vec.Mat4.invert(pvInvM, pvInvM);

    var startCoord = goog.vec.Vec3.createFloat32FromValues(x, y, 0);
    var endCoord = goog.vec.Vec3.createFloat32FromValues(x, y, 1);

    goog.vec.Mat4.multVec3Projective(pvInvM, startCoord, startCoord);
    goog.vec.Mat4.multVec3Projective(pvInvM, endCoord, endCoord);

    ray = goog.vec.Vec3.createFloat32();
    goog.vec.Vec3.subtract(endCoord, startCoord, ray);

    var t = (worldZ - startCoord[2]) / ray[2];
    goog.vec.Vec3.scale(ray, t, ray);
    goog.vec.Vec3.add(startCoord, ray, ray);
    return ray;
};

Spinteny.prototype.setDragHandler = function() {
    var thisObj = this;
    this.drag.startHandler = function(event) { thisObj.dragStart(event); };
    this.drag.moveHandler = function(event) { thisObj.dragMove(event); };
    this.drag.endHandler = function(event) { thisObj.dragEnd(event); };

    this.canvas.addEventListener("mousedown", this.drag.startHandler);
    this.canvas.addEventListener("touchstart", this.drag.startHandler);
};

Spinteny.prototype.dragStart = function(event) {
    if (event.touches) {
        this.drag.start = { x: event.touches[0].clientX,
                            y: event.touches[0].clientY };
    } else {
        this.drag.start = {x: event.clientX, y: event.clientY};
    }
    
    this.drag.org =
        Math.floor( ( ( this.drag.start.y 
                        - goog.style.getClientPosition(this.canvas).y )
                      / this.containerSize.height )
                    * (this.genomeCount) ); //TODO: actually do this right

    this.drag.org = Math.min(this.drag.org, this.genomeCount - 1);
    this.drag.initRotation = this.rotations[this.drag.org];

    this.canvas.addEventListener("mousemove", this.drag.moveHandler);
    this.canvas.addEventListener("mouseout", this.drag.endHandler);
    this.canvas.addEventListener("mouseup", this.drag.endHandler);
    this.canvas.addEventListener("touchmove", this.drag.moveHandler);
    this.canvas.addEventListener("touchend", this.drag.endHandler);

    this.canvas.removeEventListener("mousedown", this.drag.startHandler);
    this.canvas.removeEventListener("touchstart", this.drag.startHandler);
    event.preventDefault();
};

Spinteny.prototype.dragEnd = function(event) {
    this.canvas.removeEventListener("mousemove", this.drag.moveHandler);
    this.canvas.removeEventListener("mouseout", this.drag.endHandler);
    this.canvas.removeEventListener("mouseup", this.drag.endHandler);
    this.canvas.removeEventListener("touchmove", this.drag.moveHandler);
    this.canvas.removeEventListener("touchend", this.drag.endHandler);

    this.canvas.addEventListener("mousedown", this.drag.startHandler);
    this.canvas.addEventListener("touchstart", this.drag.startHandler);

    event.preventDefault();
};

Spinteny.prototype.dragMove = function(event) {
    var clientDeltaX;
    if (event.touches) {
        clientDeltaX = event.touches[0].clientX - this.drag.start.x;
    } else {
        clientDeltaX = event.clientX - this.drag.start.x;
    }
    var nearDeltaX =
        (clientDeltaX / this.containerSize.width)
        * this.widthAt(this.cameraDistance);
    var angle = (nearDeltaX / (this.genomeRadius * this.zoomFactor));
    // the % (2 * Math.PI) here is to keep the rotation magnitudes
    // low; otherwise the numeric precision issues are worse
    this.rotations[this.drag.org] =
        (this.drag.initRotation + angle) % (2 * Math.PI);
    this.updateTransform(this.drag.org);
    
    this.draw();

    event.preventDefault();
};

function triangleBarycenters(length) {
    var centers = new Float32Array([1, 0, 0,
                                    0, 1, 0,
                                    0, 0, 1]);
    var allCenters = new Float32Array(length);
    for (var i = 0; i < length; i += 9) {
        allCenters.set(centers, i);
    }
    return allCenters;
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
 *             vertex: vertex x, y, z coord (model space)
 *             normal: vertex normal x, y, z
 *             org: organism ID [ 0..(# of genomes) )
 *         }
 *     twists:
 *         {
 *             start: position of the start of the line this vertex is on
 *                    (model space)
 *             end: position of the end of the line this vertex is on
 *                  (model space)
 *             otherStart: position of the start of the other edge of the twist
 *                        (model space)
 *             otherEnd: position of the end of the other edge of the twist
 *                      (model space)
 *             prevOrg: organism ID of the previous genome
 *             nextOrg: organism ID of the next genome
 *             distance: how far along the line from start to end to
 *                       place this vertex
 *         }
 * }
 *
 * "anchors" visually represent the span of an LCB on a particular genome.
 * "twists" visually represent the connections between related anchors.
 * "twists" has extra vertex attributes because we need to calculate
 * vertex positions and normals within the vertex shader for the twists.
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
    // anchors are planar quadrilaterals
    // using drawArrays means 6 vertices per face
    var anchVertCount = 6 * numAnchors;
    var anchors = {
        // multiplying by 3 here because each vertex/normal has
        // (x, y, z) components
        vertex: new Float32Array(3 * anchVertCount),
        normal: new Float32Array(3 * anchVertCount),
        // ideally, org would be an unsigned int (or short or byte) array,
        // but GLSL (in webGL 1.0) doesn't allow those types for attributes
        org:    new Float32Array(anchVertCount)
    };

    var rowsPerTwist = 11;
    // two vertices per row, plus 3 vertices for degenerate triangles and
    // two to align the ends with the anchors
    var vertsPerTwist = (rowsPerTwist * 2) + 6;
    var twistVertTotal = vertsPerTwist * numTwists;
    var twists = {
        start:      new Float32Array(3 * twistVertTotal),
        end:        new Float32Array(3 * twistVertTotal),
        otherStart: new Float32Array(3 * twistVertTotal),
        otherEnd:   new Float32Array(3 * twistVertTotal),
        prevOrg:    new Float32Array(twistVertTotal),
        nextOrg:    new Float32Array(twistVertTotal),
        distance:   new Float32Array(twistVertTotal),
        normMult:   new Float32Array(twistVertTotal)
    };

    var curAnchor = 0;
    var twistVert = 0;
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

            // the top two vertices for twistCorners are the bottom two
            // from the previous anchor, and the bottom two vertices for
            // twistCorners are the top two vertices from the next anchor
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
            twistCorners = [oldAnchorVerts[2], oldAnchorVerts[3],
                            anchorVerts[0], anchorVerts[1]];

            var prevOrg = org;
            org = block[match][orgId];

            // Using GL.TRIANGLE_STRIP is a win, in that it lowers the
            // number of vertices that have to be processed by the
            // hardware, and tells the hardware more about shared
            // edges between triangles.  The downside is that if we
            // want to draw a bunch of disconnected strips with one
            // draw call, we need to insert extra vertices at the
            // beginning and end of each twist, which make degenerate
            // (zero-area) triangles that don't show up.  Otherwise
            // all the twists would be connected into one long strip.
            //
            // the first extra vertex is just a repeat of the first
            // vertex of the twist.
            this.addTwistVert(twists,
                              oldAnchorVerts[2], anchorVerts[0],
                              oldAnchorVerts[3], anchorVerts[1],
                              prevOrg, org, twistVert++,
                              0,  1);
            // these two are to make the top edge align with oldAnchor
            this.addTwistVert(twists,
                              oldAnchorVerts[2], anchorVerts[0],
                              oldAnchorVerts[3], anchorVerts[1],
                              prevOrg, org, twistVert++, 
                              0, 1);
            this.addTwistVert(twists,
                              oldAnchorVerts[3], anchorVerts[1],
                              oldAnchorVerts[2], anchorVerts[0],
                              prevOrg, org, twistVert++,
                              0, -1);
            // see dev-notes/twist-distance.png
            var triHalfHeight = 1 / (2 * rowsPerTwist);
            for (var row = 0; row < rowsPerTwist; row++) {
                this.addTwistVert(twists,
                                  oldAnchorVerts[2], anchorVerts[0],
                                  oldAnchorVerts[3], anchorVerts[1],
                                  prevOrg, org, twistVert++, 
                                  ((row * 2) + 1) * triHalfHeight,
                                  1);
                this.addTwistVert(twists,
                                  oldAnchorVerts[3], anchorVerts[1],
                                  oldAnchorVerts[2], anchorVerts[0],
                                  prevOrg, org, twistVert++,
                                  ((row * 2) + 2) * triHalfHeight,
                                  -1);
            }
            // this additional vertex is to make the bottom edge align with
            // the next anchor
            this.addTwistVert(twists,
                              oldAnchorVerts[2], anchorVerts[0],
                              oldAnchorVerts[3], anchorVerts[1],
                              prevOrg, org, twistVert++, 
                              1, -1);
            // for degenerate triangles, repeat the last vertex of the twist
            this.addTwistVert(twists,
                              oldAnchorVerts[3], anchorVerts[1],
                              oldAnchorVerts[2], anchorVerts[0],
                              prevOrg, org, twistVert++,
                              1, -1);
            this.addTwistVert(twists,
                              oldAnchorVerts[3], anchorVerts[1],
                              oldAnchorVerts[2], anchorVerts[0],
                              prevOrg, org, twistVert++,
                              1, -1);

            trianglesForQuad(anchorVerts, anchors.vertex, curAnchor * 18);
            repeat(org, anchors.org, curAnchor * 6, 6);
            curAnchor++;
        }
    }

    calcArrayFaceNormals(anchors.vertex, anchors.normal);

    return { anchors: anchors, twists: twists };
};

Spinteny.prototype.addTwistVert = function(twists, start, end,
                                           otherStart, otherEnd,
                                           prevOrg, nextOrg,
                                           vertIndex, distance, normMult) {
    twists.start.set(start, vertIndex * 3);
    twists.end.set(end, vertIndex * 3);
    twists.otherStart.set(otherStart, vertIndex * 3);
    twists.otherEnd.set(otherEnd, vertIndex * 3);
    twists.prevOrg[vertIndex] = prevOrg;
    twists.nextOrg[vertIndex] = nextOrg;
    twists.distance[vertIndex] = distance;
    twists.normMult[vertIndex] = normMult;
};

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
