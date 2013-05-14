/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

goog.provide('Spinteny');
goog.exportSymbol('Spinteny', Spinteny);

goog.require("goog.asserts");
goog.require("goog.dom");
goog.require("goog.debug.Logger");
goog.require("goog.events");
goog.require("goog.style");
goog.require("goog.ui.Component");
goog.require("goog.ui.Slider");
goog.require("goog.vec.Mat4");
goog.require("goog.vec.Vec3");

var alignVertShader = 
    [
        "uniform    mat4    orgTransforms[8];",
        "uniform    mat4    viewMatrix;",
        "uniform    mat4    projMatrix;",

        "attribute  vec3    thisStart;",
        "attribute  vec3    thisEnd;",
        "attribute  vec3    otherStart;",
        "attribute  vec3    otherEnd;",
        "attribute  float   prevOrg;",
        "attribute  float   nextOrg;",
        "attribute  float   dist;",
        "attribute  float   normMult;",
        "varying vec3 vNormalWorld;",

        "void main(void) {",
        // twist vertices are all along the left and right sides of
        // each twist.  thisStart and thisEnd are the start and end
        // points of the line this vertex is on.  otherStart and
        // otherEnd are the start and end points of the opposide side
        // of the twist.  dist ([0, 1]) is how far along those
        // lines this vertex goes.
        //
        // To calculate the normal, we take the cross product of the
        // line we're on and a vector across the twist to the point at
        // the same dist along the other edge of the twist.
        // see dev-notes/twist-shader.png

        "  mat4 startTrans = orgTransforms[int(prevOrg)];",
        "  mat4 endTrans = orgTransforms[int(nextOrg)];",
        "  vec4 thisStartWorld = startTrans * vec4(thisStart, 1.0);",
        "  vec4 thisEndWorld = endTrans * vec4(thisEnd, 1.0);",
        "  vec4 thisPosWorld = thisStartWorld + (dist * (thisEndWorld - thisStartWorld));",

        "  vec4 otherStartWorld = startTrans * vec4(otherStart, 1.0);",
        "  vec4 otherEndWorld = endTrans * vec4(otherEnd, 1.0);",
        "  vec4 otherPosWorld = otherStartWorld + (dist * (otherEndWorld - otherStartWorld));",

        "  vec4 thisVecWorld = thisEndWorld - thisStartWorld;",
        "  vec4 otherVecWorld = otherPosWorld - thisPosWorld;",
        "  vNormalWorld = normMult * cross(thisVecWorld.xyz, otherVecWorld.xyz);",

        "  gl_Position = projMatrix * viewMatrix * vec4(thisPosWorld.xyz, 1.0);",
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

        "varying  vec3   vNormalWorld;",

        "void main() {",
        "    float z = gl_FragCoord.z / gl_FragCoord.w;",
        "    float fogFactor = (fogRange.y - z) / (fogRange.y - fogRange.x);",
        "    fogFactor = clamp(fogFactor, 0.0, 1.0);",

        "    vec3 normal = normalize(vNormalWorld);",
        // flip normal for rear-facing fragments, because we want our
        // materials to be double-sided
        "    normal = normal * ( -1.0 + 2.0 * float( gl_FrontFacing ) );",

        "    float dirDiffuse = max(dot(normal, dirVector), 0.0);",
        "    vec3 dirLighted = dirWeight * dirDiffuse * dirColor;",
        "    vec3 ambLighted = ambWeight * ambColor;",
        "    vec3 color = dirLighted + ambLighted;",
        "    gl_FragColor = vec4(mix(vec3(1.0), color, fogFactor), 1.0);",
        "}"
    ].join( "\n" );

/** @constructor */
Spinteny = function(container, orgChroms, LCBs) {
    this.container = goog.dom.getElement(container);
    this.containerSize = goog.style.getSize(this.container);

    this.context = new GLOW.Context({
        width: this.containerSize.width,
        height: this.containerSize.height,
        antialias: true
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
    // see dev-notes/genome-radius.png
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

    var alignVerts = this.LCBsToVertices(LCBs);

    var alignShaderInfo = {
        vertexShader: alignVertShader,
        fragmentShader: alignFragShader,

        data: {
            // uniforms

            "orgTransforms": { value: this.orgTransformFlat },
            "viewMatrix": { value: this.viewMatrix },
            "projMatrix": { value: this.projMatrix },
            "fogRange": { value: this.fogRange },

            "ambColor": { value: this.ambColor },
            "ambWeight": { value: this.ambWeight },
            "dirColor": { value: this.dirColor },
            "dirWeight": { value: this.dirWeight },
            "dirVector": { value: this.dirVector },

            // attributes

            "thisStart": alignVerts.thisStart,
            "thisEnd": alignVerts.thisEnd,
            "otherStart": alignVerts.otherStart,
            "otherEnd": alignVerts.otherEnd,
            "prevOrg": alignVerts.prevOrg,
            "nextOrg": alignVerts.nextOrg,
            "dist": alignVerts.dist,
            "normMult": alignVerts.normMult
        },
        primitives: GL.TRIANGLE_STRIP
    };

    this.aligns = new GLOW.Shader(alignShaderInfo);

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
    this.aligns.draw();
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

/**
 * takes a single genome match in the form of an array:
 * [orgID, chrID, start, end]
 * and returns a quadrilateral that covers that region of the
 * spinteny cylinder in the given quad
 *
 * the current convention is that each org's drawing area has
 * top at 0 and bottom at 1, which will get transformed by that org's
 * transformation matrix in the vertex shader
 */
Spinteny.prototype.matchToQuad = function(match, quad) {
    this.mappers[match[0]].toSpatial(match[1], match[2], 0, quad.topLeft);
    this.mappers[match[0]].toSpatial(match[1], match[3], 0, quad.topRight);
    this.mappers[match[0]].toSpatial(match[1], match[2], 1, quad.botLeft);
    this.mappers[match[0]].toSpatial(match[1], match[3], 1, quad.botRight);
}

/**
 * takes an array of locally collinear blocks of the format:
 * [ [orgID, chrID, start, end], [orgID, chrID, start, end], ... ]
 *
 * and returns an object with
 * { 
 *     thisStart: position of the start of the line this vertex is on
 *                (model space)
 *     thisEnd: position of the end of the line this vertex is on
 *              (model space)
 *     otherStart: position of the start of the other edge of the twist
 *                 (model space)
 *     otherEnd: position of the end of the other edge of the twist
 *               (model space)
 *     prevOrg: organism ID of the previous genome
 *     nextOrg: organism ID of the next genome
 *     dist: how far along the line from start to end to place this vertex
 * }
 *
 * see dev-notes/twist-shader.png
 */
Spinteny.prototype.LCBsToVertices = function(blocks) {
    // "anchors" visually represent the span of an LCB on a particular genome.
    // "twists" visually represent the connections between related anchors.
    var numAnchors = 0;
    var numTwists = 0;
    for (var i = 0; i < blocks.length; i++) {
        // there's an anchor for each block,
        numAnchors += blocks[i].length;
        // and there's a twist between related anchors
        numTwists += blocks[i].length - 1;
    }

    var rowsPerTwist = 12;
    // two vertices per row, plus four to align the ends with the anchors
    var vertsPerTwist = (rowsPerTwist * 2) + 4;
    var alignVertTotal =
        ( (vertsPerTwist * numTwists) // verts for twists
          + (numAnchors * 4)          // verts for anchors
          + (2 * blocks.length) );    // verts for degen triangles
    var alignVerts = {
        thisStart:      new Float32Array(3 * alignVertTotal),
        thisEnd:        new Float32Array(3 * alignVertTotal),
        otherStart: new Float32Array(3 * alignVertTotal),
        otherEnd:   new Float32Array(3 * alignVertTotal),
        prevOrg:    new Float32Array(alignVertTotal),
        nextOrg:    new Float32Array(alignVertTotal),
        dist:   new Float32Array(alignVertTotal),
        normMult:   new Float32Array(alignVertTotal)
    };

    var curAnchor = 0;
    // vertIndex has an extra level of indirection because the value is
    // modified by a callee
    var vertIndex = new Int32Array([0]);
    var orgId = 0, chrId = 1, start = 2, end = 3;

    var anchorBuf = new Float32Array(12);
    var oldAnchorBuf = new Float32Array(12);
    var anchorVerts = {
        topLeft: anchorBuf.subarray(0, 3),
        topRight: anchorBuf.subarray(3, 6),
        botLeft: anchorBuf.subarray(6, 9),
        botRight: anchorBuf.subarray(9, 12)
    };
    var oldAnchorVerts = {
        topLeft: oldAnchorBuf.subarray(0, 3),
        topRight: oldAnchorBuf.subarray(3, 6),
        botLeft: oldAnchorBuf.subarray(6, 9),
        botRight: oldAnchorBuf.subarray(9, 12)
    };

    // the top two vertices for twistCorners are the bottom two
    // from the previous anchor, and the bottom two vertices for
    // twistCorners are the top two vertices from the next anchor
    //
    //       TL           TR
    //       |  oldAnchor |
    //       BL___________BR
    //       /            /
    //      /    twist   /
    //     /___________ /
    //    TL           TR
    //    |   anchor   |
    //    BL           BR
    var twistCorners = {
        topLeft: oldAnchorVerts.botLeft,
        topRight: oldAnchorVerts.botRight,
        botLeft: anchorVerts.topLeft,
        botRight: anchorVerts.topRight
    };

    for (var blockIdx = 0; blockIdx < blocks.length; blockIdx++) {
        var block = blocks[blockIdx];
        var match = 0;

        // take the genome-space match and convert it into
        // four 3d-space vertex positions for the quadrilateral
        // that will visually represent the match
        this.matchToQuad(block[match], anchorVerts);

        var org = block[match][orgId];

        // Using GL.TRIANGLE_STRIP is a win, in that it lowers the
        // number of vertices that have to be processed by the
        // hardware, and tells the hardware more about shared
        // edges between triangles.  The downside is that if we
        // want to draw a bunch of disconnected strips with one
        // draw call, we need to insert extra vertices at the
        // beginning and end of each strip, which make degenerate
        // (zero-area) triangles that don't show up.  Otherwise
        // all the strips would be connected into one long strip.

        // this first vert repeats for the degenerate triangle separating
        // this block from the previous one.
        this.addAlignVert(alignVerts,
                          anchorVerts.topLeft, anchorVerts.botLeft,
                          anchorVerts.topRight, anchorVerts.botRight,
                          org, org, vertIndex, 0, 1);

        // four verts for the anchor
        this.addAlignVertPair(alignVerts, anchorVerts, org, org,
                              vertIndex, 0, 0);
        this.addAlignVertPair(alignVerts, anchorVerts, org, org,
                              vertIndex, 1, 1);
        
        for (match = 1; match < block.length; match++) {
            oldAnchorBuf.set(anchorBuf);
            this.matchToQuad(block[match], anchorVerts);

            var prevOrg = org;
            org = block[match][orgId];

            // see dev-notes/twist-distance.png

            // two verts to align the top edge of the twist with oldAnchor
            this.addAlignVertPair(alignVerts, twistCorners, prevOrg, org,
                                  vertIndex, 0, 0);
            var triHalfHeight = 1 / (2 * rowsPerTwist);
            for (var row = 0; row < rowsPerTwist; row++) {
                this.addAlignVertPair(alignVerts, twistCorners, prevOrg, org,
                                      vertIndex,
                                      ((row * 2) + 1) * triHalfHeight,
                                      ((row * 2) + 2) * triHalfHeight);
            }
            // two verts to align the bottom edge of the twist with the anchor
            this.addAlignVertPair(alignVerts, twistCorners, prevOrg, org,
                                  vertIndex, 1, 1);
            // four verts for the anchor
            this.addAlignVertPair(alignVerts, anchorVerts, org, org,
                                  vertIndex, 0, 0);
            this.addAlignVertPair(alignVerts, anchorVerts, org, org,
                                  vertIndex, 1, 1);

        }
        // for degenerate triangles, repeat the last vertex of the block
        this.addAlignVert(alignVerts,
                          anchorVerts.topRight, anchorVerts.botRight,
                          anchorVerts.topLeft, anchorVerts.botLeft,
                          org, org, vertIndex, 1, -1);
    }
    //console.log([vertIndex, alignVerts.normMult.length]);

    return alignVerts;
};

Spinteny.prototype.addAlignVertPair = function(alignVerts, corners, prevOrg, 
                                               nextOrg, vertIndex, 
                                               distance1, distance2) {
    this.addAlignVert(alignVerts,
                      corners.topLeft, corners.botLeft,
                      corners.topRight, corners.botRight,
                      prevOrg, nextOrg, vertIndex,
                      distance1, 1);
    this.addAlignVert(alignVerts,
                      corners.topRight, corners.botRight,
                      corners.topLeft, corners.botLeft,
                      prevOrg, nextOrg, vertIndex,
                      distance2, -1);
};

Spinteny.prototype.addAlignVert = function(alignVerts,
                                           start, end,
                                           otherStart, otherEnd,
                                           prevOrg, nextOrg,
                                           vertIndex, distance, normMult) {
    var vert = vertIndex[0];
    vertIndex[0] += 1;
    alignVerts.thisStart.set(start, vert * 3);
    alignVerts.thisEnd.set(end, vert * 3);
    alignVerts.otherStart.set(otherStart, vert * 3);
    alignVerts.otherEnd.set(otherEnd, vert * 3);
    alignVerts.prevOrg[vert] = prevOrg;
    alignVerts.nextOrg[vert] = nextOrg;
    alignVerts.dist[vert] = distance;
    alignVerts.normMult[vert] = normMult;
};

/**
 * maps between genomic and spatial coordinates on a cylinder
 * @param chroms array of {start, end, name (optional)} objects
 * @param axis {Vec3} axis vector at the center of the cylinder
 * @param origin {Vec3} vector from the axis to the zero point
 * @param padding spacing between chroms (radians) (optional)
 *
 * @constructor
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
    this.tmpMat4 = goog.vec.Mat4.createFloat32();
    this.tmpVec3 = goog.vec.Vec3.createFloat32();

    var partialSum = 0;
    for (var i = 0; i < chroms.length; i++) {
        var name = ("name" in chroms[i]) ? chroms[i]["name"] : i;
        this.byName["" + name] = i;
        this.partialSums.push(partialSum);
        partialSum += chroms[i]["end"] - chroms[i]["start"];
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
CylMapper.prototype.toSpatial = function(index, base, distance, result) {
    distance = (distance === undefined) ? 0 : distance;
    //var index = this.byName["" + chrom];
    var angle = 
        ( ( ( this.partialSums[index] + ( base - this.chroms[index]["start"] ) )
            / this.totalLength )
          * ( 2 * Math.PI ) )
        + (this.padding * index);
    // scale the angle down to leave room for padding
    angle *= ( (2 * Math.PI)
               / ( (2 * Math.PI) 
                   + ( this.padding * ( this.chroms.length ) ) ) );
    goog.vec.Mat4.makeRotate(this.tmpMat4, angle,
                             this.rotAxis[0], this.rotAxis[1], this.rotAxis[2]);
    goog.vec.Mat4.multVec3NoTranslate(this.tmpMat4, this.origin, result);

    goog.vec.Vec3.scale(this.axis, distance, this.tmpVec3);
    goog.vec.Vec3.add(result, this.tmpVec3, result);

    return result;
};

//TODO: CylMapper.prototype.fromSpatial = function(spatialPos){}

