for f in output/*.pbrt; do 
    echo "Rendering $f";
    pbrt $f; 
done
