#include <SDL/SDL.h>
#define NO_SDL_GLEXT
#include <GL/gl.h>
#include <assert.h>

#define SCRW (640)
#define SCRH (480)
#define SCRBPP (32)

#define SDL_GL_SET_ATTR(name, value) assert(!SDL_GL_SetAttribute(SDL_GL_ ## name, value))

#define DECLARE_GL_PROC(type, name) type name = 0
/*#define DEFINE_GL_PROC(type, name) assert((name = (type)SDL_GL_GetProcAddress(#name)))*/
#define DEFINE_GL_PROC(type, name) \
	do { \
		union { type f; void *p; } u; \
		assert((u.p = SDL_GL_GetProcAddress(#name))); \
		name = u.f; \
	} while (0)

#define DO_GL_PROCS(WHAT) \
	WHAT(PFNGLCREATEPROGRAMOBJECTARBPROC, glCreateProgramObjectARB); \
	WHAT(PFNGLDELETEOBJECTARBPROC, glDeleteObjectARB); \
	WHAT(PFNGLCREATESHADEROBJECTARBPROC, glCreateShaderObjectARB); \
	WHAT(PFNGLSHADERSOURCEARBPROC, glShaderSourceARB); \
	WHAT(PFNGLCOMPILESHADERARBPROC, glCompileShaderARB); \
	WHAT(PFNGLGETOBJECTPARAMETERIVARBPROC, glGetObjectParameterivARB); \
	WHAT(PFNGLATTACHOBJECTARBPROC, glAttachObjectARB); \
	WHAT(PFNGLGETINFOLOGARBPROC, glGetInfoLogARB); \
	WHAT(PFNGLLINKPROGRAMARBPROC, glLinkProgramARB); \
	WHAT(PFNGLUSEPROGRAMOBJECTARBPROC, glUseProgramObjectARB); \
	WHAT(PFNGLGETUNIFORMLOCATIONARBPROC, glGetUniformLocationARB); \
	WHAT(PFNGLUNIFORM1FARBPROC, glUniform1fARB); \
	WHAT(PFNGLUNIFORM1IARBPROC, glUniform1iARB);

DO_GL_PROCS(DECLARE_GL_PROC)
/*
const GLcharARB *vertex_shader_code =
"void main() { gl_FrontColor = gl_Color; gl_Position = ftransform(); }"
;

const GLcharARB *fragment_shader_code =
"void main() { gl_FragColor = gl_Color; }"
;
*/
size_t
getFileSize(FILE *fp)
{
	size_t s;
	fseek(fp, 0, SEEK_END);
	s = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	return s;
}

const GLcharARB *
loadShader(const char *fname)
{
	FILE *fp = fopen(fname, "rb");
	size_t s;
	GLcharARB *b;
	assert(fp);
	s = getFileSize(fp);
	b = (GLcharARB *)malloc(s + 1);
	assert(fread(b, 1, s, fp) == s);
	b[s] = '\0';
	fclose(fp);
	return b;
}

int
main()
{
	GLenum shader_prog, shader_vert, shader_frag;
	int i;
	char buf[1024];
	const GLcharARB *vertex_shader_code, *fragment_shader_code;

	assert((vertex_shader_code = loadShader("shader.vert")));
	assert((fragment_shader_code = loadShader("shader.frag")));
	assert(!SDL_Init(SDL_INIT_EVERYTHING));
	atexit(SDL_Quit);
	SDL_GL_SET_ATTR(RED_SIZE, 8);
	SDL_GL_SET_ATTR(GREEN_SIZE, 8);
	SDL_GL_SET_ATTR(BLUE_SIZE, 8);
	SDL_GL_SET_ATTR(DEPTH_SIZE, 16);
	SDL_GL_SET_ATTR(DOUBLEBUFFER, 1);
	assert(SDL_SetVideoMode(SCRW, SCRH, SCRBPP, SDL_OPENGL /*| SDL_FULLSCREEN*/));
/*	puts((const char *)glGetString(GL_EXTENSIONS));*/
	DO_GL_PROCS(DEFINE_GL_PROC)
	assert((shader_prog = glCreateProgramObjectARB()));
	assert((shader_vert = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB)));
	assert((shader_frag = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB)));
	i = strlen(vertex_shader_code);
	printf("vert len: %d\n", i);
	glShaderSourceARB(shader_vert, 1, &vertex_shader_code, &i);
	i = strlen(fragment_shader_code);
	printf("frag len: %d\n", i);
	glShaderSourceARB(shader_frag, 1, &fragment_shader_code, &i);
	glCompileShaderARB(shader_vert);
	glCompileShaderARB(shader_frag);
	glAttachObjectARB(shader_prog, shader_vert);
	glAttachObjectARB(shader_prog, shader_frag);
	glLinkProgramARB(shader_prog);
	i = 0; glGetInfoLogARB(shader_vert, sizeof(buf) - 1, &i, buf); buf[i] = 0; printf("vert error: \"%s\"\n", buf);
	i = 0; glGetInfoLogARB(shader_frag, sizeof(buf) - 1, &i, buf); buf[i] = 0; printf("frag error: \"%s\"\n", buf);
	i = 0; glGetInfoLogARB(shader_prog, sizeof(buf) - 1, &i, buf); buf[i] = 0; printf("prog error: \"%s\"\n", buf);
	glUseProgramObjectARB(shader_prog);
	while (!SDL_GetKeyState(0)[SDLK_ESCAPE]) {
		SDL_PumpEvents();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glBegin(GL_QUADS);
		glVertex3f(-1.0f, 1.0f, 0.0f);
		glVertex3f(1.0f, 1.0f, 0.0f);
		glVertex3f(1.0f, -1.0f, 0.0f);
		glVertex3f(-1.0f, -1.0f, 0.0f);
		glEnd();
		SDL_GL_SwapBuffers();
	}
	return 0;
}
