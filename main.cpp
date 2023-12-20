namespace LRay {
// A class declaration representing 3d object
class _3DObject;

// A class declaration representing camera
class Camera;

// A class declaration to represent vector of 4 dimensions
class Vec4;

// A ray
class Ray;
}; // namespace LRay

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