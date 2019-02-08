# writevdb <input> <output> <sample_size> <voxel_size>
#rm -r output/*
#writevdb models/wStraight.10.pbrt output/wStraight.10.voxelsize1.vdb 10 1.0
#writevdb models/wStraight.10.pbrt output/wStraight.10.voxelsize0_1.vdb 10000 0.1
#writevdb models/wStraight.10.pbrt output/wStraight.10.voxelsize2.vdb 10 2.0
#writevdb models/wStraight.10.pbrt output/wStraight.10.voxelsize0_2.vdb 10000 0.2
writevdb models/wStraight.pbrt output/wStraight.voxelsize1.vdb 10 1.0
writevdb models/wCurly.pbrt output/wCurly.voxelsize1.vdb 10 1.0

