#version 330 core

in vec4 color;
out vec4 fragColor;
in vec3 v_Color;

uniform float offset;

void main(void)
{
	// Set the color of this fragment to the interpolated color
	// value computed by the rasterizer.

	//fragColor = color;
	fragColor = vec4(v_Color.r, v_Color.g, v_Color.b, 1);
}
