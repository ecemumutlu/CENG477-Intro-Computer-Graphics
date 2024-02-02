#version 330 core

uniform mat4 modelingMatrix;
uniform mat4 viewingMatrix;
uniform mat4 projectionMatrix;


layout(location=0) in vec3 inVertex;
layout(location=1) in vec3 inNormal;
layout(location=2) in vec3 inColor;
layout(location=3) in vec2 texCoord;

uniform vec3 cubeColor;

out vec4 fragWorldPos;
out vec3 fragWorldNor;
out vec3 v_Color;

void main(void)
{
	// Compute the world coordinates of the vertex and its normal.
	// These coordinates will be interpolated during the rasterization
	// stage and the fragment shader will receive the interpolated
	// coordinates.

	v_Color = inColor;
	fragWorldPos = modelingMatrix * vec4(inVertex, 1);
	fragWorldNor = inverse(transpose(mat3x3(modelingMatrix))) * inNormal;

    gl_Position = projectionMatrix * viewingMatrix * modelingMatrix * vec4(inVertex, 1);
}

