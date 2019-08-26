rm output/*.pbrt


# Environment map
prepare --input woman-head.scene --output output/deon_04samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=deon,pixelsamples=4,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/deon_08samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=deon,pixelsamples=8,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/deon_16samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=deon,pixelsamples=16,useDistantLight=#,useAreaLight=#,useInfiniteLight=
#prepare --input woman-head.scene --output output/deon_32samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=deon,pixelsamples=32,useDistantLight=#,useAreaLight=#,useInfiniteLight=

prepare --input woman-head.scene --output output/uniform_04samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=4,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/uniform_08samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=8,useDistantLight=#,useAreaLight=#,useInfiniteLight=
prepare --input woman-head.scene --output output/uniform_16samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=16,useDistantLight=#,useAreaLight=#,useInfiniteLight=
#prepare --input woman-head.scene --output output/uniform_256samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=256,useDistantLight=#,useAreaLight=#,useInfiniteLight=
#prepare --input woman-head.scene --output output/uniform_1024samples_venice.pbrt --propertyfiles base.properties moving_around_y.properties hdr-settings/venice.props --properties image_output_basename=,material_name=brown_hair,model_hair_filename=wCurly,shader=dualscattering,samplingMethod=uniform,pixelsamples=1024,useDistantLight=#,useAreaLight=#,useInfiniteLight=





#--- BEFORE ----

#prepare --input woman-head.scene --output output/brunette_distantlight_around_x.pbrt --propertyfiles base.properties lighting_around_x.properties --properties image_output_basename=brunette_distantlight_around_x,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#,useInfiniteLight=#

#showing caustic width
#prepare --input woman-head.scene --output output/brunette_distantlight_causticwidth.pbrt --propertyfiles base.properties causticwidth.properties --properties image_output_basename=brunette_glints,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#


#prepare --input woman-head.scene --output output/ref_brunette_distantlight_around_x.pbrt --propertyfiles base.properties distantlight_around_x.properties --properties image_output_basename=ref_brunette_distantlight_around_x,material_name=brown_hair_ref,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#
#prepare --input woman-head.scene --output output/blonde_distantlight_around_y.pbrt --propertyfiles base_production.properties distantlight_around_y.properties --properties material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#
#prepare --input woman-head.scene --output output/brunette_distantlight_around_y_elevated.pbrt --propertyfiles base.properties distantlight_around_y_elevated.properties --properties image_output_basename=brunette_distantlight_around_y_elevated,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#


# ---------------------------------------------------
# BRUNETTE MODEL (STRAIGHT HAIR)
# ---------------------------------------------------
# distant light
#prepare --input woman-head.scene --output output/brunette_distantlight_around_y.pbrt --propertyfiles base.properties lighting_around_y.properties --properties image_output_basename=brunette_distantlight_around_y,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#,useInfiniteLight=#
#prepare --input woman-head.scene --output output/brunette_distantlight_around_y_elevated.pbrt --propertyfiles base.properties lighting_around_y_elevated.properties --properties image_output_basename=brunette_distantlight_around_y_elevated,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#,useInfiniteLight=#
#prepare --input woman-head.scene --output output/brunette_distantlight_around_x.pbrt --propertyfiles base.properties lighting_around_x.properties --properties image_output_basename=brunette_distantlight_around_x,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#,useInfiniteLight=#

# area light
#prepare --input woman-head.scene --output output/brunette_arealight_around_y.pbrt --propertyfiles base.properties lighting_around_y.properties --properties image_output_basename=brunette_arealight_around_y,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=#,useAreaLight=,useInfiniteLight=#
#prepare --input woman-head.scene --output output/brunette_arealight_around_y_elevated.pbrt --propertyfiles base.properties lighting_around_y_elevated.properties --properties image_output_basename=brunette_arealight_around_y_elevated,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=#,useAreaLight=,useInfiniteLight=
#prepare --input woman-head.scene --output output/brunette_arealight_around_x.pbrt --propertyfiles base.properties lighting_around_x.properties --properties image_output_basename=brunette_arealight_around_x,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=#,useAreaLight=,useInfiniteLight=#

#
# moving viewer around
#
######prepare --input woman-head.scene --output output/ds_brunette_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ms_brunette_move_around_y,material_name=brown_hair,model_hair_filename=wStraight,shader=dualscattering,useDistantLight=,useAreaLight=#,useInfiniteLight=#



#prepare --input woman-head.scene --output output/ds_brunette_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ds_brunette_move_around_y,material_name=brown_hair,model_hair_filename=wStraight,shader=dualscattering,useDistantLight=,useAreaLight=#,useInfiniteLight=#,hairRadius=0.1

#prepare --input woman-head.scene --output output/ds_hardlyabsorbing_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ds_hardlyabsorbing_move_around_y,material_name=hardly_absorbing_hair,model_hair_filename=wStraight,shader=dualscattering,useDistantLight=,useAreaLight=#,useInfiniteLight=#

#prepare --input woman-head.scene --output output/ds_littleabsorbing_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ds_littleabsorbing_move_around_y,material_name=little_absorbing_hair,model_hair_filename=wStraight,shader=dualscattering,useDistantLight=,useAreaLight=#,useInfiniteLight=#

#prepare --input woman-head.scene --output output/ds_absorbing_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=ds_absorbing_move_around_y,material_name=absorbing_hair,model_hair_filename=wStraight,shader=dualscattering,useDistantLight=,useAreaLight=#,useInfiniteLight=#

#
# rotating model
#
#prepare --input woman-head.scene --output output/brunette_rotate_around_y.pbrt --propertyfiles base.properties rotating_model_around_y.properties --properties image_output_basename=brunette_rotating_model_around_y,material_name=brown_hair,model_hair_filename=wStraight,useDistantLight=,useAreaLight=#,useInfiniteLight=
#prepare --input woman-head.scene --output output/blonde_rotate_around_y.pbrt --propertyfiles base.properties rotating_model_around_y.properties --properties image_output_basename=blonde_rotate_model_around_y,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=

# ---------------------------------------------------
# BLONDE MODEL (WAVY HAIR)
# ---------------------------------------------------
# distant light
#prepare --input woman-head.scene --output output/blonde_distantlight_around_y.pbrt --propertyfiles base.properties distantlight_around_y.properties --properties image_output_basename=blonde_distantlight_around_y,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=#
#prepare --input woman-head.scene --output output/blonde_distantlight_around_y_elevated.pbrt --propertyfiles base.properties distantlight_around_y_elevated.48.properties --properties image_output_basename=blonde_distantlight_around_y_elevated,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=#
#prepare --input woman-head.scene --output output/blonde_distantlight_around_x.pbrt --propertyfiles base.properties distantlight_around_x.properties --properties image_output_basename=blonde_distantlight_around_x,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=#

# area light
#prepare --input woman-head.scene --output output/blonde_arealight_around_y.pbrt --propertyfiles base.properties lighting_around_y.properties --properties image_output_basename=blonde_arealight_around_y,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=#,useAreaLight=,useInfiniteLight=#
#prepare --input woman-head.scene --output output/blonde_arealight_around_y_elevated.pbrt --propertyfiles base.properties lighting_around_y_elevated.properties --properties image_output_basename=blonde_arealight_around_y_elevated,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=#,useAreaLight=,useInfiniteLight=
#prepare --input woman-head.scene --output output/blonde_arealight_around_x.pbrt --propertyfiles base.properties lighting_around_x.properties --properties image_output_basename=blonde_arealight_around_x,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=#,useAreaLight=,useInfiniteLight=#

#
# moving viewer around
#
#prepare --input woman-head.scene --output output/blonde_move_around_y.pbrt --propertyfiles base.properties moving_around_y.properties --properties image_output_basename=blonde_move_around_y,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=

#
# rotating model
#
#prepare --input woman-head.scene --output output/blonde_rotate_around_y.pbrt --propertyfiles base.properties rotating_model_around_y.properties --properties image_output_basename=blonde_rotating_model_around_y,material_name=blonde_hair,model_hair_filename=wWavy,useDistantLight=,useAreaLight=#,useInfiniteLight=
