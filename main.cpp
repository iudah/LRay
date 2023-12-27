namespace LRay {
// A class declaration representing 3d object
class _3DObject;

// A class declaration representing camera
class Camera;

// A class declaration to represent vector of 4 dimensions
class Vec4;

// A class declaration to represent vector of 3 dimensions
class Vec3;

// A ray
class Ray;
}; // namespace LRay
// Vec4
class LRay::Vec4 {
protected:
  float x;
  float y;
  float z;
  float w;

public:
LRay::Vec4(float x, float y, float z, float w=0);
LRay::Vec4 &add(LRay::Vec4 &b);
LRay::Vec4 &minus(LRay::Vec4 &b);
};

class LRay::Vec3 :public LRay::Vec4{

};

class LRay::_3DObject {};

class LRay::Ray {};

class LRay::Camera {
private:
  LRay::Vec4 position;
  const int width;
  const int height;

public:
  Camera(int width = 720, int height = 960);

  // Render image
  int render();
};

int main() {
  LRay::Camera camera_0 = LRay::Camera();
  camera_0.render();
}

