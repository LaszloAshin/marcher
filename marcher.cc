#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>

class Vector {
	float x_, y_, z_, w_;

public:
	Vector(float x = 0.0f, float y = 0.0f, float z = 0.0f, float w = 1.0f)
	: x_(x)
	, y_(y)
	, z_(z)
	, w_(w)
	{}

	float x() const { return x_; }
	float y() const { return y_; }
	float z() const { return z_; }
	float w() const { return w_; }

	Vector &x(float x) { x_ = x; return *this; }
	Vector &y(float y) { y_ = y; return *this; }
	Vector &z(float z) { z_ = z; return *this; }
	Vector &w(float w) { w_ = w; return *this; }

	Vector h() const { return Vector(x_ / w_, y_ / w_, z_ / w_); }

	Vector operator*(float s) const { return Vector(x_ * s, y_ * s, z_ * s, w_); }
	Vector operator+(const Vector &v) const { return Vector(x_ / w_ + v.x_ / v.w_, y_ / w_ + v.y_ / v.w_, z_ / w_ + v.z_ / v.w_); }
	Vector operator-(const Vector &v) const { return Vector(x_ / w_ - v.x_ / v.w_, y_ / w_ - v.y_ / v.w_, z_ / w_ - v.z_ / v.w_); }

	float l() const { float x(x_ / w_), y(y_ / w_), z(z_ / w_); return sqrtf(x * x + y * y + z * z); }
	Vector norm() const { float len(l()); return Vector(x_ / len, y_ / len, z_ / len); }

	float dot(const Vector &v) const { return x_ * v.x_ + y_ * v.y_ + z_ * v.z_; }
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
	float r_, g_, b_, a_;

public:
	Pixel(float r = 0.0f, float g = 0.0f, float b = 0.0f, float a = 1.0f)
	: r_(r)
	, g_(g)
	, b_(b)
	, a_(a)
	{}

	float r() const { return r_; }
	float g() const { return g_; }
	float b() const { return b_; }
	float a() const { return a_; }

	Pixel &r(float r) { r_ = r; return *this; }
	Pixel &g(float g) { g_ = g; return *this; }
	Pixel &b(float b) { b_ = b; return *this; }
	Pixel &a(float a) { a_ = a; return *this; }

	Pixel &min(const Pixel &p);
	Pixel &max(const Pixel &p);

	Pixel operator*(float s) const { return Pixel(r_ * s, g_ * s, b_ * s, a_ * s); }
};

Pixel &
Pixel::min(const Pixel &p)
{
	if (p.r_ < r_) r_ = p.r_;
	if (p.g_ < g_) g_ = p.g_;
	if (p.b_ < b_) b_ = p.b_;
	if (p.a_ < a_) a_ = p.a_;
	return *this;
}

Pixel &
Pixel::max(const Pixel &p)
{
	if (p.r_ > r_) r_ = p.r_;
	if (p.g_ > g_) g_ = p.g_;
	if (p.b_ > b_) b_ = p.b_;
	if (p.a_ > a_) a_ = p.a_;
	return *this;
}

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
	std::pair<Pixel, Pixel> minmax() const;
};

Pixmap &
Pixmap::clear(const Pixel &p)
{
	for (size_t i = 0; i < pm_.size(); ++i) pm_[i] = p;
	return *this;
}

std::pair<Pixel, Pixel>
Pixmap::minmax()
const
{
	Pixel min(pm_[0]), max(pm_[0]);
	for (size_t i = 1; i < pm_.size(); ++i) {
		min.min(pm_[i]);
		max.max(pm_[i]);
	}
	return std::make_pair(min, max);
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
	Pixel max(pm.minmax().second);
	float maxval = max.r();
	if (max.g() > maxval) maxval = max.g();
	if (max.b() > maxval) maxval = max.b();
	ofs << unsigned(maxval * scale) << std::endl;
	for (size_t y = 0; y < pm.h(); ++y) {
		for (size_t x = 0; x < pm.w(); ++x) {
			Pixel p(pm.at(x, y) * scale);
			ofs << std::setw(4) << unsigned(p.r());
			ofs << std::setw(4) << unsigned(p.g());
			ofs << std::setw(4) << unsigned(p.b());
			ofs << std::setw(4) << "";
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

class FuncFace {
	static const float eps = 0.01f;

public:
	virtual ~FuncFace() {}

	virtual float f(float x, float z, Vector *n) const = 0;
	Vector dx(float x, float z) const { return Vector(2.0f * eps, f(x + eps, z, 0) - f(x - eps, z, 0), 0.0f); }
	Vector dz(float x, float z) const { return Vector(0.0f, f(x, z + eps, 0) - f(x, z - eps, 0), 2.0f * eps); }
	Vector n(float x, float z) const { return (dz(x, z) % dx(x, z)).norm(); }
};

template <class T>
class Interpolator {
public:
	virtual ~Interpolator() {}

	virtual T f(float x) const = 0;
	virtual T df(float x) const = 0;
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
		return (1.0f - x) * 6.0f * x;
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
		return (1.0f + (x - 2.0f) * x) * 30.0f * x * x;
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

class MyFuncFace : public FuncFace {
	static float noise(int x) { return 1.0f - ((x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0; }
	static float noise2d(int x, int z) { return noise(z * 1024 + x); }

public:
//	virtual float f(float x, float z) const { return sinf(x) * sinf(z); }
	float g(float x, float z, Vector *n) const;
	float h(float x, float z, Vector *n) const;
	virtual float f(float x, float z, Vector *n) const;
};

float
MyFuncFace::g(float x, float z, Vector *n)
const
{
	int ix = int(floor(x));
	int iz = int(floor(z));
	x -= ix, z -= iz;
	float m[4][4] = {
		{ noise2d(ix - 1, iz - 1), noise2d(ix, iz - 1), noise2d(ix + 1, iz - 1), noise2d(ix + 2, iz - 1) },
		{ noise2d(ix - 1, iz), noise2d(ix, iz), noise2d(ix + 1, iz), noise2d(ix + 2, iz) },
		{ noise2d(ix - 1, iz + 1), noise2d(ix, iz + 1), noise2d(ix + 1, iz + 1), noise2d(ix + 2, iz + 1) },
		{ noise2d(ix - 1, iz + 2), noise2d(ix, iz + 2), noise2d(ix + 1, iz + 2), noise2d(ix + 2, iz + 2) }
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
MyFuncFace::h(float x, float z, Vector *n)
const
{
	int ix = int(floor(x));
	int iz = int(floor(z));
	float ox = x, oz = z;
	x -= ix, z -= iz;
	float m[2][2] = {
		{ noise2d(ix, iz), noise2d(ix + 1, iz) },
		{ noise2d(ix, iz + 1), noise2d(ix + 1, iz + 1) }
	};
	float a = SmoothstepInterpolator<float>(m[0][0], m[0][1]).f(x);
	float b = SmoothstepInterpolator<float>(m[1][0], m[1][1]).f(x);
	SmoothstepInterpolator<float> ciz(a, b);
	if (n) {
/*		const float eps = 0.01f;
		Vector dx(2.0f * eps, h(ox + eps, oz, 0) - h(ox - eps, oz, 0), 0.0f);
		Vector dz(0.0f, h(ox, oz + eps, 0) - h(ox, oz - eps, 0), 2.0f * eps);
		*n = Vector(dz % dx).norm();*/
		float c = SmootherstepInterpolator<float>(m[0][0], m[1][0]).f(z);
		float d = SmootherstepInterpolator<float>(m[0][1], m[1][1]).f(z);
		SmootherstepInterpolator<float> cix(c, d);
		*n = Vector(cix.df(x), 1.0f, ciz.df(z));
	}
	return ciz.f(z);
}

float
MyFuncFace::f(float x, float z, Vector *n)
const
{
	float y = 0.0f;
	float a = 1.0f;
	Vector no;
	for (int i = 0; i < 6; ++i) {
		float yo = a * h(x, z, n);
//		if (n) yo /= n->dot(*n);
		y += yo;
		if (n) no = no + *n;
		a *= 0.5f;
		x *= 2.0f;
		z *= 2.0f;
	}
	if (n) {
		no.y(1.0f);
		*n = no.norm();
	}
	return y;
}

class Scene {
	const FuncFace &gnd_;
	Vector light_;
	Vector cam_;

public:
	Scene(FuncFace &gnd)
	: gnd_(gnd)
	{}

	Scene &light(const Vector &light) { light_ = light; return *this; }
	Scene &cam(const Vector &cam) { cam_ = cam; return *this; }

	const Vector &cam() const { return cam_; }

	Pixel cast(const Ray &r);
};

Pixel
Scene::cast(const Ray &r)
{
	const float step = 0.02f;
	float dt = step;
	float ldt = dt;
	const float tmin = 1.0f;
	const float tmax = 10.0f;
	float ld = 0.0f;
	for (float t = tmin; t < tmax; t += dt) {
		Vector p(r.o() + r.d() * t);
		if (p.y() > 2.0f && r.d().y() > 0.0f) break;
		float h = gnd_.f(p.x(), p.z(), 0);
		float cd = h - p.y();
		if (cd > 0.0f) {
			t -= ldt * cd / (ld + cd);
			p = r.o() + r.d() * t;
			gnd_.f(p.x(), p.z(), 0);
			Vector n(gnd_.n(p.x(), p.z()));
			Vector l((light_ - p).norm());
			Vector v((cam_ - p).norm());
			Vector k(0.1f, l.dot(n));
			if (k.y() > 0.0f) {
				k.z(powf(((n * k.y() - l) * 2.0f + l).dot(v), 8.0f));
			} else {
				k.y(0.0f);
			}
//			return Pixel(h * h);
			return Pixel(1.0f, 1.0f, 1.0f) * k.dot(Vector(0.1f, 0.6f, 0.3f));
			if (n.x() < 0.0f) n.x(0.0f);
			if (n.y() < 0.0f) n.y(0.0f);
			if (n.z() < 0.0f) n.z(0.0f);
			return Pixel(n.x(), n.y(), n.z());
		}
		ldt = dt;
		ld = -cd;
		dt = t * step;
	}
	return Pixel();
}

class Tracer {
	float fov_;

public:
	Tracer(float fov) : fov_(fov) {}

	Pixmap &render(Pixmap &pm, Scene &s);
};

Pixmap &
Tracer::render(Pixmap &pm, Scene &s)
{
	Vector d;
	d.z(-1.0f / tanf((0.5f * fov_) * (M_PI / 180.0f)));
	for (int y = 0; y < pm.h(); ++y) {
		d.y(1.0f - 2.0f * y / (pm.h() - 1));
		d.y(d.y() * pm.h() / pm.w());
		for (int x = 0; x < pm.w(); ++x) {
			d.x(2.0f * x / (pm.w() - 1) - 1.0f);
			pm.at(x, y) = s.cast(Ray(s.cam(), d));
		}
		std::cout << y << std::endl;
	}
	return pm;
}

int
main()
{
	Pixmap pm(640, 480);
	MyFuncFace ff1;
	Scene s(ff1);
	s.light(Vector(0.0f, 4.0f, 0.0f));
	s.cam(Vector(0.0f, 1.0f, 0.0f));
	Netbpm("x.ppm").save(Tracer(90.0f).render(pm, s));
	return 0;
}
