#version 430

in VertexOut {
    vec3 color;
} fsIn;

layout(location=0) out vec4 colorOut;

void main() {
    colorOut = vec4(1.0f, 1.0f, 0.0f, 1.0f);
}
