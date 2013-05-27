goog.provide('CylMapper');

goog.require("goog.vec.Mat4");
goog.require("goog.vec.Vec3");

/**
 * maps between genomic and spatial coordinates on a cylinder
 * @param chroms array of {start, end, name (optional)} objects
 * @param axis {Vec3} axis vector at the center of the cylinder
 * @param origin {Vec3} vector from the axis to the zero point
 * @param padding spacing between chroms (radians) (optional)
 *
 * @constructor
 */
CylMapper = function(chroms, axis, origin, padding) {
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

