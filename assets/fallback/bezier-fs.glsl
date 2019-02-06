#version 120

uniform vec3 color;

//layout(location=0) out vec4 colorOut;

void main() {
    gl_FragColor = vec4(color, 1.0f);
    //colorOut = vec4(color, 1.0f);
}
