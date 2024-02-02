#version 330 core

// All of the following variables could be defined in the OpenGL
// program and passed to this shader as uniform variables. This
// would be necessary if their values could change during runtim.
// However, we will not change them and therefore we define them 
// here for simplicity.

vec3 I = vec3(1, 1, 1);          // point light intensity
vec3 Iamb = vec3(0.8, 0.8, 0.8); // ambient light intensity
vec3 kd = vec3(0.8, 0.8, 0.8);     // diffuse reflectance coefficient
vec3 ka = vec3(0.3, 0.3, 0.3);   // ambient reflectance coefficient
vec3 ks = vec3(0.8, 0.8, 0.8);   // specular reflectance coefficient
vec3 lightPos = vec3(5, 5, 5);   // light position in world coordinates

uniform vec3 eyePos;
uniform vec3 color;
uniform float offset;

in vec4 fragWorldPos;
in vec3 fragWorldNor;
in vec3 v_Color;

out vec4 fragColor;

vec3 calculateColor(vec3 pos, float scale, float offset){
	bool x = int((pos.x + 70.0f ) * scale) % 2 != 0;
    bool y = int((pos.y + 70.0f ) * scale) % 2 != 0;
    bool z = int((pos.z + offset) * scale) % 2 != 0;
	bool xorXY = x != y;

	if (xorXY != z)
		return vec3(0.0f, 0.0f, 1.0f); //blue
	else
		return vec3(1.0f, 1.0f, 1.0f); //white
}


void main(void)
{
	// Compute lighting. We assume lightPos and eyePos are in world
	// coordinates. fragWorldPos and fragWorldNor are the interpolated
	// coordinates by the rasterizer.
	v_Color;
	vec3 new_Color = v_Color;

	new_Color = calculateColor(fragWorldPos.xyz, 0.115f, offset);
	vec3 L = normalize(lightPos - vec3(fragWorldPos));
	vec3 V = normalize(eyePos - vec3(fragWorldPos));
	vec3 H = normalize(L + V);
	vec3 N = normalize(fragWorldNor);

	float NdotL = dot(N, L); // for diffuse component
	float NdotH = dot(N, H); // for specular component

	vec3 diffuseColor = I  * new_Color * kd * max(0, NdotL);
	vec3 specularColor = I * ks * pow(max(0, NdotH), 100) ;
	vec3 ambientColor = Iamb * ka * new_Color ;

	fragColor = vec4(diffuseColor + specularColor + ambientColor, 1);
}
