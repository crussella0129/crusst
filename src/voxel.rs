//! Voxelization: sample an SDF onto a regular 3D grid.
//!
//! The [`VoxelGrid`] stores signed-distance values (as `f32`) at each voxel
//! center. Construction is parallelized with rayon for large grids.

use nalgebra::Vector3;

/// A regular 3D grid of signed-distance values.
pub struct VoxelGrid {
    /// Number of voxels along each axis `[nx, ny, nz]`.
    pub resolution: [usize; 3],
    /// Edge length of each cubic voxel.
    pub voxel_size: f64,
    /// World-space position of the grid's minimum corner (voxel `[0,0,0]`
    /// has its center at `origin + voxel_size * 0.5`).
    pub origin: Vector3<f64>,
    /// Flat array of SDF values in row-major order `[x][y][z]`.
    /// Index = `ix * ny * nz + iy * nz + iz`.
    pub data: Vec<f32>,
}

impl VoxelGrid {
    /// Convert a world-space point to the flat index of the nearest voxel.
    ///
    /// Points outside the grid are clamped to the boundary.
    pub fn index_at(&self, world: Vector3<f64>) -> usize {
        let [nx, ny, nz] = self.resolution;
        let rel = world - self.origin;
        let ix = ((rel.x / self.voxel_size) as isize).clamp(0, nx as isize - 1) as usize;
        let iy = ((rel.y / self.voxel_size) as isize).clamp(0, ny as isize - 1) as usize;
        let iz = ((rel.z / self.voxel_size) as isize).clamp(0, nz as isize - 1) as usize;
        ix * ny * nz + iy * nz + iz
    }

    /// Return the world-space center of the voxel at grid indices `(ix, iy, iz)`.
    pub fn world_at(&self, ix: usize, iy: usize, iz: usize) -> Vector3<f64> {
        self.origin
            + Vector3::new(
                (ix as f64 + 0.5) * self.voxel_size,
                (iy as f64 + 0.5) * self.voxel_size,
                (iz as f64 + 0.5) * self.voxel_size,
            )
    }
}
