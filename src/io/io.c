#include <stdio.h>
#include <zot.h>

bool parse_objs(FILE *f)
{
    char ln[256];
    while (fgets(ln, sizeof(ln), f))
    {
    }
    return true;
}

bool load_objs(const char *fpath)
{
    FILE *f = fopen(fpath, "r");

    if (!f)
    {
        LOG("File `%s` could not be opened", fpath);
        return false;
    }

    return parse_objs(f);
}
bool load_camera_motion(const char *fpath)
{
    return true;
}
bool load_camera_keyframes(const char *fpath) { return true; }