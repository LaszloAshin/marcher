varying vec4 pos;
uniform float time;
const float eps = 1e-6;
vec3 L1;
vec4 CO;
float t;
vec3 C, V, N;

void
sphere(in vec4 col, in vec3 d, in float r)
{
	float D = pow(dot(d, V), 2) - dot(V, V) * (dot(d, d) - r * r);
	if (D > eps) {
		D = sqrt(D);
		float t1 = (-dot(d, V) + D) / dot(V, V);
		float t2 = (-dot(d, V) - D) / dot(V, V);
		if (t2 < eps) t2 = t1;
		if (t2 > eps && t2 < t) {
			t = t2;
			N = d + V * t2;
			CO = col;
		}
	}
}

void
plane(in vec4 col, in vec3 q, in float d)
{
	float a = dot(V, q);
	if (abs(a) > eps) {
		float t2 = (d - dot(C, q)) / a;
		if (t2 > eps && t2 < t) {
			t = t2;
			N = q;
			CO = col;
		}
	}
}

void
scene()
{
		sphere(vec4(0.0, 1.0, 0.0, 0.0), C - vec3(2.0*cos(time), sin(time), 4.0+sin(time*0.87)), 1.0);
		sphere(vec4(0.0, 0.0, 1.0, 0.0), C - vec3(-2.0*cos(time), 2.0*sin(time * 0.95), 3.0+sin(time*0.77)), 1.0);
		plane(vec4(1.0, 0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0), -1.0);
//XXX:	plane(vec4(1.0, 1.0, 0.0, 0.0), vec3(0.0, -1.0, 0.0), 1.0);
}

void
main()
{
	L1 = vec3(50.0 * sin(0.1 * time), 10.0, 50.0 * cos(0.1 * time));
	C = vec3(0.0, 0.0, 0.0);
	V = vec3(pos.x, pos.y, 1.0);
	vec4 a = vec4(0.0, 0.0, 0.0, 1.0);
	for (int i = 0; i < 4; ++i) {
		t = 1000.0;
		scene();
		if (t > 999.0) {
			break;
		}
		N = normalize(N);
		vec3 p = C + V * t;
		vec3 L = normalize(L1 - p);
		vec3 Cn = normalize(C - p);
		vec3 R = 2.0 * (N * dot(Cn, N)) - Cn;
		vec3 m = vec3(0.3, 0.6, 0.1);
		vec3 k = vec3(max(0.0, dot(L, N)), pow(max(0.0, dot(L, R)), 32), 1.0);
		vec4 col = CO;
		C = p + L * 1e-3;
		V = L;
		t = 1000.0;
		scene();
		if (t < 999.0) {
			// in shadow
			k.x = k.y = 0.0;
		}
		a += col * dot(m, k) * (1.0 / pow(2.0, i));
		C = p + R * 1e-3;
		V = R;
	}
	gl_FragColor = clamp(a, vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0));
}
