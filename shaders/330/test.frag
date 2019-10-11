#version 330 core
out vec4 FragColor;
in vec3 sharedPosition;

uniform int isBox;
uniform vec3 myLightPosition;

void main()
{
	if (isBox == 0) {
		vec3 N;
		vec4 ambient = vec4(49, 100, 183, 255) / 255;
		vec4 diffuse_light = vec4(171, 202, 252, 255) / 255;
		vec3 lightDir = normalize(sharedPosition - myLightPosition);

		N.xy = gl_PointCoord * 2.0 - vec2(1.0); 

		float mag = dot(N.xy, N.xy);

		if (mag > 1.0) discard;   // kill pixels outside circle
		N.z = sqrt(1.0-mag);

		// calculate lighting
		float diffuse = max(0.0, dot(lightDir, N));


		FragColor = normalize(diffuse_light * diffuse + ambient);
	} else {
		FragColor = vec4(1.0, 1.0, 1.0, 1.0);
	}
} 