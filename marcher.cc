/**
 * Raymarcher by Laszlo Ashin 2o11.
 * g++ -c marcher.cc -march=native -O3 -fomit-frame-pointer -pipe -Wall -pedantic -pthread
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>

// http://rrrola.wz.cz/inv_sqrt.html
float
invsqrt(float x)
{
	union { float f; unsigned u; } y = { x };
	y.u = 0x5F1FFFF9ul - (y.u >> 1);
	return 0.703952253f * y.f * (2.38924456f - x * y.f * y.f);
}

class Vector {
	float x_, y_, z_;

public:
	explicit
	Vector(float x = 0.0f, float y = 0.0f, float z = 0.0f)
	: x_(x)
	, y_(y)
	, z_(z)
	{}

	float x() const { return x_; }
	float y() const { return y_; }
	float z() const { return z_; }

	Vector &x(float x) { x_ = x; return *this; }
	Vector &y(float y) { y_ = y; return *this; }
	Vector &z(float z) { z_ = z; return *this; }

	Vector operator*(float s) const { return Vector(x_ * s, y_ * s, z_ * s); }
	Vector operator+(const Vector &v) const { return Vector(x_ + v.x_, y_ + v.y_, z_ + v.z_); }
	Vector operator-(const Vector &v) const { return Vector(x_ - v.x_, y_ - v.y_, z_ - v.z_); }

	float dot(const Vector &v) const { return x_ * v.x_ + y_ * v.y_ + z_ * v.z_; }

	float il() const { return invsqrt(dot(*this)); }
	Vector norm() const { return *this * il(); }
	float len() const { return 1.0f / il(); }

	Vector operator%(const Vector &v) const;
};

Vector
Vector::operator%(const Vector &v)
const
{
	return Vector(
		y_ * v.z_ - z_ * v.y_,
		z_ * v.x_ - x_ * v.z_,
		x_ * v.y_ - y_ * v.x_
	);
}

class Pixel {
	float r_, g_, b_;

public:
	explicit
	Pixel(float r = 0.0f, float g = 0.0f, float b = 0.0f)
	: r_(r)
	, g_(g)
	, b_(b)
	{}

	float r() const { return r_; }
	float g() const { return g_; }
	float b() const { return b_; }

	Pixel &r(float r) { r_ = r; return *this; }
	Pixel &g(float g) { g_ = g; return *this; }
	Pixel &b(float b) { b_ = b; return *this; }

	Pixel operator*(float s) const { return Pixel(r_ * s, g_ * s, b_ * s); }
	Pixel operator+(const Pixel &p) const { return Pixel(r_ + p.r_, g_ + p.g_, b_ + p.b_); }
	Pixel operator-(const Pixel &p) const { return Pixel(r_ - p.r_, g_ - p.g_, b_ - p.b_); }
	float dot(const Pixel &p) const { return r_ * p.r_ + g_ * p.g_ + b_ * p.b_; }

	float d(const Pixel &p) const { Pixel d(*this - p); return d.dot(d); }
};

class Pixmap {
	size_t w_, h_;
	std::vector<Pixel> pm_;

	size_t o(size_t x, size_t y) const { return y * w_ + x; }

public:
	Pixmap(size_t w, size_t h, const Pixel &p = Pixel())
	: w_(w)
	, h_(h)
	, pm_(w_ * h_, p)
	{}

	size_t w() const { return w_; }
	size_t h() const { return h_; }

	Pixel &at(int x, int y) { return pm_[o(x, y)]; }
	const Pixel &at(int x, int y) const { return pm_[o(x, y)]; }

	Pixmap &clear(const Pixel &p = Pixel());
};

Pixmap &
Pixmap::clear(const Pixel &p)
{
	for (size_t i = 0; i < pm_.size(); ++i) pm_[i] = p;
	return *this;
}

class Netbpm {
	const std::string &fn_;

public:
	Netbpm(const std::string &fn)
	: fn_(fn)
	{}

	Netbpm &save(const Pixmap &pm, float scale = 255.0f);
};

Netbpm &
Netbpm::save(const Pixmap &pm, float scale)
{
	std::ofstream ofs(fn_.c_str());
	ofs << "P3" << std::endl;
	ofs << pm.w() << " " << pm.h() << std::endl;
	ofs << int(scale) << std::endl;
	for (size_t y = 0; y < pm.h(); ++y) {
		for (size_t x = 0; x < pm.w(); ++x) {
			Pixel p(pm.at(x, y) * scale);
			ofs << int(p.r()) << " " << int(p.g()) << " " << int(p.b()) << " ";
		}
		ofs << std::endl;
	}
	return *this;
}

class Ray {
	Vector o_, d_;

public:
	Ray(const Vector &o, const Vector &d)
	: o_(o)
	, d_(d.norm())
	{}

	const Vector &o() const { return o_; }
	const Vector &d() const { return d_; }
};

float
clamp(float x, float l = 0.0f, float h = 1.0f)
{
	if (x < l) {
		return l;
	} else if (x > h) {
		return h;
	} else {
		return x;
	}
}

template <class T>
class Interpolator {
public:
	virtual ~Interpolator() {}

	virtual T f(float x) const = 0;
	virtual T df(float x) const = 0;
};

template <class T>
class LinearInterpolator : public Interpolator<T> {
	const T p0_, p1_;

public:
	LinearInterpolator(const T &p0, const T &p1)
	: p0_(p0)
	, p1_(p1)
	{}

	virtual T f(float x) const { return p0_ * (1.0f - x) + p1_ * x; }
	virtual T df(float x) const { 0.0f; }
};

template <class T>
class SmoothstepInterpolator : public Interpolator<T> {
	const T p0_, p1_;

public:
	SmoothstepInterpolator(const T &p0, const T &p1)
	: p0_(p0)
	, p1_(p1)
	{}

	virtual T
	f(float x)
	const
	{
		float t = (3.0f - 2.0f * x ) * x * x;
		return p0_ * (1.0f - t) + p1_ * t;
	}

	virtual T
	df(float x)
	const
	{
		float t = (1.0f - x) * 6.0f * x;
		return p0_ * (1.0f - t) + p1_ * t; // XXX: ez szar
	}
};

template <class T>
class SmootherstepInterpolator : public Interpolator<T> {
	const T p0_, p1_;

public:
	SmootherstepInterpolator(const T &p0, const T &p1)
	: p0_(p0)
	, p1_(p1)
	{}

	virtual T
	f(float x)
	const
	{
		float t = (10.0f + (6.0f * x - 15.0f) * x) * x * x * x;
		return p0_ * (1.0f - t) + p1_ * t;
	}

	virtual T
	df(float x)
	const
	{
		float t = (1.0f + (x - 2.0f) * x) * 30.0f * x * x;
		return p0_ * (1.0f - t) + p1_ * t; // XXX: ez szar
	}
};

template <class T>
class CubicInterpolator : public Interpolator<T> {
	const T p0_, p1_, p2_, p3_;

public:
	CubicInterpolator(const T &p0, const T &p1, const T &p2, const T &p3)
	: p0_(p0)
	, p1_(p1)
	, p2_(p2)
	, p3_(p3)
	{}

	virtual T
	f(float x)
	const
	{
		T y = (p1_ - p2_) * 1.5f + (p3_ - p0_) * 0.5f;
		y *= x;
		y += p0_ - p1_ * 2.5 + p2_ * 2.0f - p3_ * 0.5f;
		y *= x;
		y += (p2_ - p0_) * 0.5f;
		return y * x + p1_;
	}

	virtual T
	df(float x)
	const
	{
		T y = (p1_ - p2_) * 4.5f + (p3_ - p0_) * 1.5f;
		y *= x;
		y += p0_ * 2.0f - p1_ * 5.0f + p2_ * 4.0f - p3_;
		y *= x;
		return y + (p2_ - p0_) * 0.5f;
	}
};

class FuncFace {
public:
	virtual ~FuncFace() {}

	virtual float f(float x, float z) const = 0;
	virtual Pixel color(const Vector &p, const Vector &k, const Vector &n) const = 0;
};

class FuncFaceDerivator {
	const FuncFace &f_;
	const float eps_;

public:
	FuncFaceDerivator(const FuncFace &face, float eps = 0.01f)
	: f_(face)
	, eps_(eps)
	{}

	Vector dx(float x, float z) const { return Vector(2.0f * eps_, f_.f(x + eps_, z) - f_.f(x - eps_, z), 0.0f); }
	Vector dz(float x, float z) const { return Vector(0.0f, f_.f(x, z + eps_) - f_.f(x, z - eps_), 2.0f * eps_); }
	Vector n(float x, float z) const { return (dz(x, z) % dx(x, z)).norm(); }
};

class PerlinNoise {
	static float noise1f(int x) { return 1.0f - ((x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0; }
	static float noise2f(int x, int y) { return noise1f(y * 1024 + x); }
	static float noise3f(int x, int y, int z) { return noise1f((z * 1024 + y) * 1024 + x); }

public:
	static float bicubic(float x, float z, Vector *n);
	static float bilinear(float x, float z);
	static float trilinear(float x, float y, float z);
};

float
PerlinNoise::bicubic(float x, float z, Vector *n)
{
	int ix = int(floor(x));
	int iz = int(floor(z));
	x -= ix, z -= iz;
	float m[4][4] = {
		{ noise2f(ix - 1, iz - 1), noise2f(ix, iz - 1), noise2f(ix + 1, iz - 1), noise2f(ix + 2, iz - 1) },
		{ noise2f(ix - 1, iz), noise2f(ix, iz), noise2f(ix + 1, iz), noise2f(ix + 2, iz) },
		{ noise2f(ix - 1, iz + 1), noise2f(ix, iz + 1), noise2f(ix + 1, iz + 1), noise2f(ix + 2, iz + 1) },
		{ noise2f(ix - 1, iz + 2), noise2f(ix, iz + 2), noise2f(ix + 1, iz + 2), noise2f(ix + 2, iz + 2) }
	};
	float a = CubicInterpolator<float>(m[0][0], m[0][1], m[0][2], m[0][3]).f(x);
	float b = CubicInterpolator<float>(m[1][0], m[1][1], m[1][2], m[1][3]).f(x);
	float c = CubicInterpolator<float>(m[2][0], m[2][1], m[2][2], m[2][3]).f(x);
	float d = CubicInterpolator<float>(m[3][0], m[3][1], m[3][2], m[3][3]).f(x);
	CubicInterpolator<float> ciz(a, b, c, d);
	if (n) {
		float e = CubicInterpolator<float>(m[0][0], m[1][0], m[2][0], m[3][0]).f(z);
		float f = CubicInterpolator<float>(m[0][1], m[1][1], m[2][1], m[3][1]).f(z);
		float g = CubicInterpolator<float>(m[0][2], m[1][2], m[2][2], m[3][2]).f(z);
		float h = CubicInterpolator<float>(m[0][3], m[1][3], m[2][3], m[3][3]).f(z);
		CubicInterpolator<float> cix(e, f, g, h);
		*n = Vector(-cix.df(x), 1.0f, -ciz.df(z));
	}
	return ciz.f(z);
}

float
PerlinNoise::bilinear(float x, float z)
{
	int ix = int(floor(x));
	int iz = int(floor(z));
	x -= ix, z -= iz;
	float m[2][2] = {
		{ noise2f(ix, iz), noise2f(ix + 1, iz) },
		{ noise2f(ix, iz + 1), noise2f(ix + 1, iz + 1) }
	};
	float a = SmoothstepInterpolator<float>(m[0][0], m[0][1]).f(x);
	float b = SmoothstepInterpolator<float>(m[1][0], m[1][1]).f(x);
	SmoothstepInterpolator<float> ciz(a, b);
	return ciz.f(z);
}

float
PerlinNoise::trilinear(float x, float y, float z)
{
	int ix = int(floor(x));
	int iy = int(floor(y));
	int iz = int(floor(z));
	x -= ix, y -= iy, z -= iz;
	float a = SmoothstepInterpolator<float>(noise3f(ix, iy, iz), noise3f(ix + 1, iy, iz)).f(x);
	float b = SmoothstepInterpolator<float>(noise3f(ix, iy + 1, iz), noise3f(ix + 1, iy + 1, iz)).f(x);
	float c = SmoothstepInterpolator<float>(noise3f(ix, iy, iz + 1), noise3f(ix + 1, iy, iz + 1)).f(x);
	float d = SmoothstepInterpolator<float>(noise3f(ix, iy + 1, iz + 1), noise3f(ix + 1, iy + 1, iz + 1)).f(x);
	float e = SmoothstepInterpolator<float>(a, b).f(y);
	float f = SmoothstepInterpolator<float>(c, d).f(y);
	return SmoothstepInterpolator<float>(e, f).f(z);
}

class Ground : public FuncFace {
public:
	virtual float f(float x, float z) const;
	virtual Pixel color(const Vector &p, const Vector &k, const Vector &n) const;
};

float
Ground::f(float x, float z)
const
{
	float y = 0.0f;
	float a = 2.0f;
	for (int i = 0; i < 8; ++i) {
		float m = 1.0f / (1.0f + 2.0f * !i);
		y += a * PerlinNoise::bilinear(x * m, z * m);
		a *= 0.3f * m;
		x *= 2.0f;
		z *= 2.0f;
	}
	return y;
}

Pixel
Ground::color(const Vector &p, const Vector &k, const Vector &n)
const
{
	float r = 0.0f;
	float a = 1.0f;
	float x = 4.0f * p.x();
	float z = 4.0f * p.z();
	for (int i = 0; i < 6; ++i) {
		r += a * PerlinNoise::bilinear(x, z);
		a *= 0.5f;
		x *= 2.0f;
		z *= 2.0f;
	}
	float t = 1.0f - (1.0f + r * r) / 3.0f;
	Pixel sand(Pixel(0.8f, 0.6f, 0.4f) * (t * k.dot(Vector(0.3f, 0.7f, 0.0f))));
	Pixel grass(Pixel(0.2f, 0.8f, 0.2f) * (t * k.dot(Vector(0.2f, 0.7f, 0.1f))));
	Pixel snow(Pixel(1.0f, 1.0f, 1.0f) * k.dot(Vector(0.8f, 0.2f, 0.0f)));
	Pixel pix(SmoothstepInterpolator<Pixel>(sand, grass).f(pow(clamp(n.y() + 0.1f), 4.0f)));
	return SmoothstepInterpolator<Pixel>(pix, snow).f(pow(clamp(p.y() + r, 0.0f, 1.0f), 8.0f));
}

class Water : public FuncFace {
public:
	virtual float f(float x, float z) const;
	virtual Pixel color(const Vector &p, const Vector &k, const Vector &n) const;
};

float
Water::f(float x, float z)
const
{
	float y = -1.7f;
	float a = 0.01f;
	x *= 10.0f, z *= 10.0f;
	for (int i = 0; i < 3; ++i) {
		y += a * PerlinNoise::bilinear(x, z);
		a *= 0.5f;
		x *= 2.0f;
		z *= 2.0f;
	}
	return y;
}

Pixel
Water::color(const Vector &p, const Vector &k, const Vector &n)
const
{
	return Pixel(0.1f, 0.5f, 1.0f) * k.dot(Vector(0.1f, 0.4f, 0.5f));
}

class Clouds {
public:
	float f(float x, float y, float z) const;
};

float
Clouds::f(float x, float y, float z)
const
{
	float f = 0.0f;
	float a = 1.0f;
//	x /= 10.0f, y /= 10.0f, z /= 10.0f;
	for (int i = 0; i < 3; ++i) {
		f += a * PerlinNoise::trilinear(x, y, z);
		a *= 0.5f;
		x *= 2.0f;
		y *= 2.0f;
		z *= 2.0f;
	}
	return f;
}

class Scene {
	Ground gnd_;
	Water wtr_;
	Clouds cld_;
	Vector light_;
	Vector cam_;

public:
	Scene()
	: gnd_()
	, wtr_()
	, cld_()
	, light_()
	, cam_()
	{}

	Scene &light(const Vector &light) { light_ = light; return *this; }
	Scene &cam(const Vector &cam) { cam_ = cam; return *this; }

	const Vector &cam() const { return cam_; }

	bool isInShadow(const Vector &v) const;
	Pixel cast(const Ray &r) const;
};

bool
Scene::isInShadow(const Vector &v)
const
{
	Ray r(v, (light_ - v).norm());
	const float step = 0.2f;
	float dt = step;
	const float tmin = 0.0f;
	const float tmax = (light_ - v).len();
	for (float t = tmin; t < tmax; t += dt) {
		dt = (1.0f + t) * step;
		Vector p(r.o() + r.d() * t);
		if (p.y() > 2.0f) {
			if (r.d().y() > 0.0f) break;
		}
		if (gnd_.f(p.x(), p.z()) > p.y()) {
			return true;
		}
	}
	return false;
}

Pixel
Scene::cast(const Ray &r)
const
{
	const float step = 0.01f;
	float dt = step;
	float ldt = dt;
	const float tmin = 1.0f;
	const float tmax = 25.0f;
	float ldg = 0.0f;
	float ldw = 0.0f;
	float blu = r.d().y() * r.d().y();
	blu = 0.5f - 0.5f * blu;
	Pixel pix(blu, blu, 1.0f);
	float sh = 0.0f;
	float d = 0.0f;
	for (float t = tmin; t < tmax; t += dt) {
		dt = t * step; // non-linear steps
		Vector p(r.o() + r.d() * t);
		float pd = tmax / 30.0f * dt;
		if (p.y() > 2.0f) {
			pd = dt * cld_.f(p.x() - 0.99f * r.o().x(), p.y() - 0.99f * r.o().y(), p.z());
		} else if (isInShadow(p)) {
			sh += pd;
		}
		d += pd;
		if (p.y() > 2.0f) {
			if (r.d().y() > 0.0f) continue;
		}
		float hg = gnd_.f(p.x(), p.z());
		float hw = wtr_.f(p.x(), p.z());
		float cdg = hg - p.y();
		float cdw = hw - p.y();
		if (fmaxf(cdg, cdw) > 0.0f) {
			if (cdg > cdw) {
				t -= ldt * cdg / (ldg + cdg);
			} else {
				t -= ldt * cdw / (ldw + cdw);
			}
			p = r.o() + r.d() * t;
			const FuncFace *ff;
			if (wtr_.f(p.x(), p.z()) > gnd_.f(p.x(), p.z())) {
				ff = &wtr_;
			} else {
				ff = &gnd_;
			}
			bool inShad = isInShadow(p - (r.d() * ldt));
			Vector n(FuncFaceDerivator(*ff).n(p.x(), p.z()));
			Vector l(((inShad ? cam_ : light_) - p).norm());
			Vector v((cam_ - p).norm());
			Vector k(1.0f, l.dot(n));
			if (k.y() > 0.0f) {
				k.z(powf(((n * k.y() - l) * 2.0f + l).dot(v), 8.0f));
			} else {
				k.y(0.0f);
			}
			if (inShad) k = k * 0.5f;
//			return Pixel(h * h);
//			float f = fmaxf(t / tmax, 0.75f * (-1.2f - p.y()));
//			float f = t / tmax;
//			return SmoothstepInterpolator<Pixel>(ff->color(k), Pixel(0.5f, 0.5f, 0.5f)).f(powf(clamp(f, 0.0f, 1.0f), 0.5f));
			pix = ff->color(p, k, n);
			break;
/*			if (n.x() < 0.0f) n.x(0.0f);
			if (n.y() < 0.0f) n.y(0.0f);
			if (n.z() < 0.0f) n.z(0.0f);
			return Pixel(n.x(), n.y(), n.z());*/
		}
		ldt = dt;
		ldg = -cdg;
		ldw = -cdw;
	}
	float col = 1.0f - clamp(sh / d, 0.0f, 1.0f);
//	float col = 1.0f;
	return SmoothstepInterpolator<Pixel>(pix, Pixel(col, col, col)).f(pow(clamp(0.25f * d / tmax, 0.0f, 1.0f), 0.25f));
}

class Tracer {
	float fov_;

public:
	Tracer(float fov) : fov_(fov) {}

	static float halton(int base, int n);
	Ray genRay(const Pixmap &pm, float x, float y, const Scene &s) const;
	Pixmap &render(Pixmap &pm, Scene &s);
};

float
Tracer::halton(int base, int n)
{
	float ret = 0.0f;
	int b = base;

	while (n) {
		ret += (float)(n % base) / b;
		b *= base;
		n /= base;
	}
	return ret;
}

Ray
Tracer::genRay(const Pixmap &pm, float x, float y, const Scene &s)
const
{
	return Ray(
		s.cam(),
		Vector(
			2.0f * x / pm.w() - 1.0f,
			(1.0f - 2.0f * y / pm.h()) * pm.h() / pm.w(),
			-1.0f / tanf((0.5f * fov_) * (M_PI / 180.0f))
		)
	);
}

Pixmap &
Tracer::render(Pixmap &pm, Scene &s)
{
	const float EDGE_LIMIT = 0.01f;
	const int OVERSAMPLE = 8;

	std::vector<Pixel> pl(pm.w() + 1);
	for (size_t x = 0; x <= pm.w(); ++x) {
		pl[x] = s.cast(genRay(pm, x, 0.0f, s));
	}
	for (size_t y = 0; y < pm.h(); ++y) {
		Pixel pp(s.cast(genRay(pm, 0.0f, 1.0f + y, s)));
		for (size_t x = 0; x < pm.w(); ++x) {
			Pixel p(s.cast(genRay(pm, 1.0f + x, 1.0f + y, s)));
			Pixel a = pl[x] + pl[x + 1] + pp + p;
			int na = 4;
			if (pl[x].d(p) > EDGE_LIMIT || pl[x + 1].d(pp) > EDGE_LIMIT) {
				for (int i = 1; i <= OVERSAMPLE; ++i) {
					a = a + s.cast(genRay(pm, halton(2, i) + x, halton(3, i) + y, s));
					++na;
				}
			}
			pm.at(x, y) = a * (1.0f / na);
			pl[x] = pp;
			pp = p;
		}
		pl[pm.w()] = pp;
		std::cout << y << std::endl;
	}
	return pm;
}

int
main()
{
	Pixmap pm(1920, 1080);
	Scene s;
	s.cam(Vector(-22.0f, 1.0f, 1.0f));
	s.light(s.cam() + Vector(-17.0f, 4.0f, -11.0f));
	Netbpm("x.ppm").save(Tracer(120.0f).render(pm, s));
	return 0;
}
