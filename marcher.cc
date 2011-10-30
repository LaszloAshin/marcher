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

	virtual float f(float x, float z) const = 0;
	Vector dx(float x, float z) const { return Vector(2.0f * eps, f(x + eps, z) - f(x - eps, z), 0.0f); }
	Vector dz(float x, float z) const { return Vector(0.0f, f(x, z + eps) - f(x, z - eps), 2.0f * eps); }
	Vector n(float x, float z) const { return (dz(x, z) % dx(x, z)).norm(); }
};

template <class T>
class Interpolator {
public:
	virtual ~Interpolator() {}

	virtual T at(float x) const = 0;
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
	at(float x)
	const
	{
		float fx = -0.5f * p0_ + 1.5f * p1_ - 1.5f * p2_ + 0.5f * p3_;
		fx *= x;
		fx += p0_ - 2.5 * p1_ + 2.0f * p2_ - 0.5f * p3_;
		fx *= x;
		fx += -0.5f * p0_ + 0.5f * p2_;
		return x * fx + p1_;
	}
};

class MyFuncFace : public FuncFace {
	static float noise(int x) { return 1.0f - ((x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0; }
	static float noise2d(int x, int z) { return noise(z * 1024 + x); }

public:
//	virtual float f(float x, float z) const { return sinf(x) * sinf(z); }
	virtual float f(float x, float z) const;
};

float
MyFuncFace::f(float x, float z)
const
{
	int ix = int(floor(x));
	int iz = int(floor(z));
	x -= ix, z -= iz;
	float a = CubicInterpolator<float>(noise2d(ix - 1, iz - 1), noise2d(ix, iz - 1), noise2d(ix + 1, iz - 1), noise2d(ix + 2, iz - 1)).at(x);
	float b = CubicInterpolator<float>(noise2d(ix - 1, iz), noise2d(ix, iz), noise2d(ix + 1, iz), noise2d(ix + 2, iz)).at(x);
	float c = CubicInterpolator<float>(noise2d(ix - 1, iz + 1), noise2d(ix, iz + 1), noise2d(ix + 1, iz + 1), noise2d(ix + 2, iz + 1)).at(x);
	float d = CubicInterpolator<float>(noise2d(ix - 1, iz + 2), noise2d(ix, iz + 2), noise2d(ix + 1, iz + 2), noise2d(ix + 2, iz + 2)).at(x);
	return CubicInterpolator<float>(a, b, c, d).at(z);
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
	float tinc = 0.01f;
	float tmin = 1.0f;
	float tmax = 50.0f;
	for (float t = tmin; t < tmax; t += tinc) {
		Vector p(r.o() + r.d() * t);
		float h = gnd_.f(p.x(), p.z());
		if (p.y() < h) {
			if (tinc > 0.1f) {
				t -= tinc;
				tinc = 0.01f;
				continue;
			}
			Vector n(gnd_.n(p.x(), p.z()));
			Vector l((light_ - p).norm());
			Vector v((cam_ - p).norm());
			Vector k(0.1f, l.dot(n));
			if (k.y() > 0.0f) {
				k.z(powf(((n * k.y() - l) * 2.0f + l).dot(v), 8.0f));
			} else {
				k.y(0.0f);
			}
			return Pixel(1.0f, 1.0f, 1.0f) * k.dot(Vector(0.1f, 0.6f, 0.3f));
		}
		tinc *= 1.1f;
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
		for (int x = 0; x < pm.w(); ++x) {
			d.x(2.0f * x / (pm.w() - 1) - 1.0f);
			pm.at(x, y) = s.cast(Ray(s.cam(), d));
		}
	}
	return pm;
}

int
main()
{
	Pixmap pm(640, 480);
	MyFuncFace ff1;
	Scene s(ff1);
	s.light(Vector(-4.0f, 4.0f, -10.0f));
	s.cam(Vector(0.0f, 4.0f, 0.0f));
	Netbpm("x.ppm").save(Tracer(90.0f).render(pm, s));
	return 0;
}
