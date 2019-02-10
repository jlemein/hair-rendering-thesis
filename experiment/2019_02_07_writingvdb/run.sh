# writevdb <input> <output> <sample_size> <voxel_size>
rm -r output/*
writevdb models/wStraight.pbrt output/wStraight.voxelsize0_5.vdb 10 0.5
#writevdb models/wStraight.pbrt output/wStraight.voxelsize5.vdb 1000 5.0
#writevdb models/wStraight.10.pbrt output/wStraight.10.voxelsize2.vdb 10 2.0
#writevdb models/wStraight.pbrt output/wStraight.voxelsize3.vdb 1000 3.0
#writevdb models/wStraight.pbrt output/wStraight.voxelsize1.vdb 10 1.0
#writevdb models/wCurly.pbrt output/wCurly.voxelsize1.vdb 10 1.0

