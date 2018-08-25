rm output/*.pbrt

# straight hair
prepare --input woman-head.scene --output output/distant_front_over_back.pbrt --propertyfiles base.properties distant_front_over_back.properties --properties material_name=brown_hair,model_hair_filename=wWavy
prepare --input woman-head.scene --output output/distant_from_back_around_CW.pbrt --propertyfiles base.properties distant_from_back_around_CW.properties --properties material_name=blonde_hair,model_hair_filename=wWavy


