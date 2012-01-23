varying vec4 pos;
const float eps = 1e-6;
vec3 L1 = vec3(10.0, 10.0, 0.0);
vec4 CO;

void
sphere(inout float t, inout vec3 N, in vec3 V, in vec4 col, in vec3 d, in float r)
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
plane(inout float t, inout vec3 N, in vec3 V, in vec4 col, in vec3 C, in vec3 q, in float d)
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

vec4
trace(in vec3 C, in vec3 V)
{
	vec4 c = vec4(0.0, 0.0, 0.0, 1.0);
	for (int i = 0; i < 4; ++i) {
		float t = 1000.0f;
		vec3 N;
		sphere(t, N, V, vec4(0.0, 1.0, 0.0, 0.0), C - vec3(2.0, -1.0, 4.0), 1.0);
		sphere(t, N, V, vec4(0.0, 0.0, 1.0, 0.0), C - vec3(-2.0, -1.0, 3.0), 1.0);
		plane(t, N, V, vec4(1.0, 0.0, 0.0, 0.0), C, vec3(0.0, 1.0, 0.0), -1.0);
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
		c += CO * dot(m, k) * (1.0 / pow(2.0, i));
		C = p + R * eps;
		V = R;
	}
	return clamp(c, vec4(0.0, 0.0, 0.0, 1.0), vec4(1.0, 1.0, 1.0, 1.0));
}

void
main()
{
	gl_FragColor = trace(vec3(0.0, 0.0, 0.0), vec3(pos.x, pos.y, 1.0));
}
