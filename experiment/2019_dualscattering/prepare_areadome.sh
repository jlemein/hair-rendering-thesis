rm output/*.pbrt

prepare --input woman-head-arealight-dome.scene --output output/ds_arealight.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ds_arealight,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,useDistantLight=#,useAreaLight=,useInfiniteLight=#
