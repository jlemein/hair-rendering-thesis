#version 120

uniform mat4 model;
uniform mat4 view;
uniform mat4 perspective;
uniform vec3 color;

//layout (location=0) in vec3 position;
varying vec3 position;

void main() {
    mat4 MVP = perspective * view * model;
    gl_Position = MVP * vec4(position, 1.0);
}
