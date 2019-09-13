rm output/*.pbrt

### ===============================================================================
### Variation of df, db from 0.0 to 1.0 (in 10 increments) lit from side
### ===============================================================================

prepare --input woman-head.scene --output output/dbdf_brunette_venice.pbrt --propertyfiles base.properties dbdf-variation.props hdr-settings/venice.props --properties image_output_basename=brown_dbdf,material_name=brown_hair,model_hair_filename=wCurly.z,shader=dualscattering,samplingMethod=deon,pixelsamples=16,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/dbdf_blonde_venice.pbrt --propertyfiles base.properties dbdf-variation.props hdr-settings/venice.props --properties image_output_basename=blonde_dbdf,material_name=blonde_hair,model_hair_filename=wCurly.z,shader=dualscattering,samplingMethod=deon,pixelsamples=16,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/dbdf_black_venice.pbrt --propertyfiles base.properties dbdf-variation.props hdr-settings/venice.props --properties image_output_basename=black_dbdf,material_name=black_hair,model_hair_filename=wCurly.z,shader=dualscattering,samplingMethod=deon,pixelsamples=16,useDistantLight=#,useAreaLight=#,useInfiniteLight=
#prepare --input woman-head.scene --output output/dbdf_uniform_venice.pbrt --propertyfiles base.properties dbdf-variation.props hdr-settings/venice.props --properties image_output_basename=dbdf,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=1,useDistantLight=#,useAreaLight=#,useInfiniteLight=
