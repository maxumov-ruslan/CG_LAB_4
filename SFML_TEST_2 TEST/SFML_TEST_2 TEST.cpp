
#define _USE_MATH_DEFINES

#include <vector>
#include <exception>
#include <cstdlib>
#include <SFML/Graphics.hpp>
#include <cmath>

using namespace sf;
using namespace std;

class Matrix {
	struct MatrixSize {
		int w;
		int h;
	} size;
	vector<double> vec;
public:
	Matrix(int h, int w, const double* coefficients) {
		size.w = w;
		size.h = h;
		vec.assign(coefficients, coefficients + w * h);
	}

	Matrix(int h, int w, initializer_list<double> coefficients) {
		size.w = w;
		size.h = h;
		vec.assign(coefficients.begin(), coefficients.end());
	}

	static Matrix I(int s, double a = 1.0) {
		double* c = new double[s * s];
		for (int y = 0; y < s; y++) {
			for (int x = 0; x < s; x++) {
				c[y * s + x] = (x == y ? a : 0.0);
			}
		}
		Matrix i(s, s, c);
		delete[] c;
		return i;
	}

	static Matrix v3m4x1(Vector3f vec) {
		double m[] = { vec.x,vec.y,vec.z,1 };
		return Matrix(4, 1, m);
	}

	Vector3f transform(Vector3f vec) {
		Matrix columnVec = Matrix::v3m4x1(vec);
		Matrix transformed = *this * columnVec;
		transformed = I(4, 1.0 / transformed.get(3, 0)) * transformed;
		return transformed.toVec();
	}

	Vector2f project(Vector3f vec) {
		Matrix columnVec = Matrix::v3m4x1(vec);
		Matrix projected = (*this) * columnVec;
		projected = I(4, 1.0 / projected.get(3, 0)) * projected;
		return projected.toVec2();
	}

	double get(int i, int j) const {
		return vec[size.w * i + j];
	}

	Vector2f toVec2() const {
		return Vector2f{ (float)get(0,0),(float)get(1,0) };
	}

	Vector3f toVec()const {
		return Vector3f((float)get(0, 0), (float)get(1, 0), (float)get(2, 0));
	}

	Matrix operator*(const Matrix& b) const {
		const Matrix& a = *this;
		if (a.size.w != b.size.h)
			throw runtime_error("Can't multiply matricies with incompatible sizes!");
		int h = a.size.h;
		int w = b.size.w;
		vector<double> v;
		for (int i = 0; i < h; i++)
			for (int j = 0; j < w; j++) {
				double sum = 0;
				for (int k = 0; k < a.size.w; k++)
					sum += a.get(i, k) * b.get(k, j);
				v.push_back(sum);
			}
		return Matrix(h, w, v.data());
	}

	Matrix operator+(const Matrix& other)const {
		if (other.size.w != size.w || other.size.h != size.h) {
			throw runtime_error("Can't add matricies with incompatible sizes!");
		}
		vector<double> c(size.w * size.h);
		for (int i = 0; i < size.h; i++)
			for (int j = 0; j < size.w; j++)
				c[i * size.w + j] = get(i, j) + other.get(i, j);
		return Matrix(size.h, size.w, c.data());
	}

	static Matrix translate(double x, double y, double z) {
		return Matrix(4, 4, { 1, 0, 0, x,
							  0, 1, 0, y,
							  0, 0, 1, z,
							  0, 0, 0, 1 });
	}

	static Matrix rotateX(double phi) {
		return Matrix(4, 4, { 1,    0,        0, 0,
							  0, cos(phi),-sin(phi),0,
							  0, sin(phi), cos(phi) ,0,
							  0, 0, 0, 1 });
	}

	static Matrix rotateY(double phi) {
		return Matrix(4, 4, { cos(phi), 0, sin(phi), 0,
								 0,     1,    0,     0,
							 -sin(phi), 0, cos(phi), 0,
								 0,     0,    0,     1 });
	}

	static Matrix rotateZ(double phi) {
		return Matrix(4, 4, { cos(phi),-sin(phi), 0, 0,
							  sin(phi), cos(phi), 0, 0,
								 0       , 0    , 1, 0,
								 0       , 0    , 0, 1 });
	}

	static Matrix scale(double x, double y, double z) {
		return Matrix(4, 4, { x, 0, 0, 0,
							  0, y, 0, 0,
							  0, 0, z, 0,
							  0, 0, 0, 1 });
	}
};
Vector2f centerCoordinates(Vector2f vec, const RenderWindow& win) {
	auto size = win.getSize();
	return Vector2f{ vec.x + size.x / 2,-vec.y + size.y / 2 };
}
int mod(int i, int n) {
	return (n + (i % n)) % n;
}
struct Shape3d {
	vector<Vector3f> vertices;
	vector<pair<int, int>> edges;
	Shape3d& transform(Matrix transformation) {
		vector<Vector3f> nVertices;
		for (const auto& v : vertices) {
			nVertices.push_back(transformation.transform(v));
		}
		vertices = nVertices;
		return *this;
	}

	void draw(RenderWindow& win, Matrix projectionMatrix = Matrix::I(4)) {
		VertexArray array;
		array.setPrimitiveType(Lines);
		for (const auto& e : edges) {
			Vector2f a = centerCoordinates(projectionMatrix.project(vertices[e.first]), win);
			Vector2f b = centerCoordinates(projectionMatrix.project(vertices[e.second]), win);
			array.append(Vertex(a, Color::White));
			array.append(Vertex(b, Color::White));
		}
		win.draw(array);
	}

	static Shape3d cube() {
		return Shape3d{
			{
				{-1,-1,-1},
				{-1,-1, 1},
				{-1, 1,-1},
				{-1, 1, 1},
				{1, -1,-1},
				{1, -1, 1},
				{1,  1,-1},
				{1,  1, 1}
			},
			{
				{0, 1},
				{0, 2},
				{0, 4},
				{1, 3},
				{1, 5},
				{2, 3},
				{2, 6},
				{3, 7},
				{4, 5},
				{4, 6},
				{5, 7},
				{7, 6}
			}
		};
	}

	static Shape3d prism(int n, double a, double h) {
		vector<Vector3f> V(2 * n);
		vector<pair<int, int>> E(3 * n);
		double r = a / sin(M_PI / n);
		for (int i = 0; i < n; i++) {
			V[2 * i] = (Vector3f(r * cos(i * 2 * M_PI / n), r * sin(i * 2 * M_PI / n), -h / 2));
			V[2 * i + 1] = (Vector3f(r * cos(i * 2 * M_PI / n), r * sin(i * 2 * M_PI / n), h / 2));
			E[3 * i] = make_pair(2 * i, 2 * i + 1);
			E[3 * i + 1] = make_pair(mod(2 * i + 2, 2 * n), 2 * i);
			E[3 * i + 2] = make_pair(mod(2 * i + 3, 2 * n), mod(2 * i + 1, 2 * n));
		}
		return Shape3d{ V, E };
	}
};
const Matrix ISOMETRIC = Matrix(4, 4, { sqrt(3),    0,   -sqrt(3),  0,
										   1,       2,      1,      0,
										sqrt(2),-sqrt(2), sqrt(2),  0,
										   0,       0,       0,  sqrt(6) });
const double d = 400.0;
const Matrix PERSPECTIVE = Matrix(4, 4, { 1, 0,   0,   0,
										  0, 1,   0,   0,
										  0, 0,   1,   d,
										  0, 0, 1 / d, 1 });
int main() {
	vector<Shape3d> shapes
	{ Shape3d::cube().transform(Matrix::scale(50, 50, 50)) ,
	Shape3d::prism(6, 50, 250).transform(Matrix::rotateX(M_PI / 2)).transform(Matrix::translate(250, 0, 0)),
	Shape3d::prism(5, 70, 200).transform(Matrix::translate(-250, 0, 0)) };

	int activeShape = 0;
	Matrix projection = PERSPECTIVE;

	RenderWindow window(VideoMode(1080, 1080), "Lab 4");

	while (window.isOpen()) {
		window.clear();
		Event event;

		while (window.pollEvent(event)) {
			if (event.type == Event::Closed) window.close();

			if (event.type == Event::KeyPressed) {
				switch (event.key.code)
				{
				case Keyboard::Down:
					shapes[activeShape].transform(Matrix::rotateZ(M_PI / 32.0));
					break;

				case Keyboard::Up:
					shapes[activeShape].transform(Matrix::rotateZ(-M_PI / 32.0));
					break;

				case Keyboard::Left:
					shapes[activeShape].transform(Matrix::rotateY(M_PI / 32.0));
					break;

				case Keyboard::Right:
					shapes[activeShape].transform(Matrix::rotateY(-M_PI / 32.0));
					break;

				case Keyboard::W:
					shapes[activeShape].transform(Matrix::translate(0, 10, 0));
					break;

				case Keyboard::S:
					shapes[activeShape].transform(Matrix::translate(0, -10, 0));
					break;

				case Keyboard::A:
					shapes[activeShape].transform(Matrix::translate(-10, 0, 0));
					break;

				case Keyboard::D:
					shapes[activeShape].transform(Matrix::translate(10, 0, 0));
					break;

				case Keyboard::Q:
					shapes[activeShape].transform(Matrix::translate(0, 0, -10));
					break;

				case Keyboard::E:
					shapes[activeShape].transform(Matrix::translate(0, 0, 10));
					break;

				case Keyboard::I:
					projection = ISOMETRIC;
					break;

				case Keyboard::P:
					projection = PERSPECTIVE;
					break;

				case Keyboard::Space:
					activeShape = (activeShape + 1) % shapes.size();

				default:
					break;
				}
			}
		}
		for (auto& shape : shapes) {
			shape.draw(window, projection);
		}
		window.display();
	}
}