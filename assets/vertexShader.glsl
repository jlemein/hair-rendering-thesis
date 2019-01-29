#version 430

uniform mat4 model;
uniform mat4 view;
uniform mat4 perspective;
uniform vec3 offsets[200];

layout (location=0) in vec3 position;
//layout (location=1) in vec3 color;

out VertexOut {
    vec3 color;
} vsOut;

void main() {
    //vsOut.color = color;
    mat4 MVP = perspective * view * model;
    gl_Position = MVP * vec4(position, 1.0);
}
