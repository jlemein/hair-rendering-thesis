#version 430

uniform vec3 color;

in VertexOut {
    vec3 vcolor;
} fsIn;

layout(location=0) out vec4 colorOut;

void main() {
    colorOut = vec4(color, 1.0f);
}
