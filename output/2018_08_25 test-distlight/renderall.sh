for f in straight-hair-with-head/*.pbrt; do 
    echo "Rendering $f";
    pbrt $f; 
done
