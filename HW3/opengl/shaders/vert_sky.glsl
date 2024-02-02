#version 330 core

layout(location=0) in vec3 inVertex;
layout(location=1) in vec2 texCoord;

out vec2 v_TexCoord;


void main()
{
   gl_Position =  vec4(inVertex.xy, 0.99, 1);
   v_TexCoord = texCoord;
};



