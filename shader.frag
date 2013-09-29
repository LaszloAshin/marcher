varying vec4 pos;
uniform float time;
const float eps = 1e-6;
vec3 V = vec3(pos.x, pos.y, 1.0);
vec3 C = vec3(0.0, 0.0, time);

float rand(float x, float y) {
    return fract(sin(12.9898 * x + 78.233 * y) * 43758.5453);
}

float noise(vec2 p) {
	vec2 i = floor(p);
	p -= i;
	p = smoothstep(0.0, 1.0, p);
	return mix(
		mix(rand(i.x, i.y), rand(i.x + 1.0, i.y), p.x),
		mix(rand(i.x, i.y + 1.0), rand(i.x + 1.0, i.y + 1.0), p.x),
		p.y
	);
}

float perlin(vec2 p) {
	float z = 0.0;
	float a = 0.5;
	for (int i = 0; i < 1; ++i) {
		z += a * noise(p);
		a *= 0.5;
		p *= 2.0;
	}
	return z;
}

void
main()
{
	float t = 1.0;
	float ody = 0.0;
	float ot = t;
	while (t < 5.0) {
		vec3 p = C + V * t;
		float f = 2.0 * perlin(p.xz) - 1.0;
		float dy = f - p.y;
		if (dy > 0.0) {
			t = (ot * dy - t * ody) / (dy - ody);
			gl_FragColor = vec4(1.0 - (t - 1.0) / 5.0, 0.0, 0.0, 1.0);
			return;
		}
		ody = dy;
		ot = t;
		t += 0.05 * t;
	}
	gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}
