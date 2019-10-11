#version 330 core
layout (location = 0) in vec3 pos;

uniform mat4 myModelviewMatrix;
uniform mat4 myProjectionMatrix;

uniform float farPlane;
uniform float maxPointSize;

uniform int isBox;

out vec3 sharedPosition;

void main() 
{
	sharedPosition = pos;
	vec4 finalPos =  myProjectionMatrix * myModelviewMatrix * vec4(pos.x, pos.y, pos.z, 1.0);

	gl_PointSize = farPlane * maxPointSize / finalPos.w;

	gl_Position = finalPos;
}

