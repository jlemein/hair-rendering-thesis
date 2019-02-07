#version 430

layout (lines_adjency) in;
layout (line_strip, max_vertices = 4) out;

// From vertex shader
//in VS_Out {
//    vec2 texCoord;
//    vec3 normal;
//    vec4 lightVector;
//    vec3 eyeVector;
//} gsIn[];


// Geometry Out
//out GS_Out {
    //vec2 texCoord;
    //vec3 normal;
    //vec4 lightVector;
    vec3 eyeVector;
//} gsOut;

void main() {
    //gsOut.lightVector = gsIn[0].lightVector;
    //gsOut.eyeVector = gsIn[0].eyeVector;

    //gsOut.texCoord = gsIn[0].texCoord;
    //gsOut.normal = gsIn[0].normal;
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    gl_Position = gl_in[2].gl_Position;
    EmitVertex();

    gl_Position = gl_in[3].gl_Position;
    EmitVertex();

    EndPrimitive();
}