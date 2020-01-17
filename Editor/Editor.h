#include <cstdint>
#include <imgui\imgui.h>
#include "imgui_impl_sdl_gl3.h"

#include <GL\glew.h>  
#include <SDL\SDL.h>
#include <chrono>
#undef main // undef to remove sdl_main

#pragma once
class VolumeData;
class Editor {
public:
	Editor(uint32_t width, uint32_t height);
	~Editor();
	bool Initialize();
	void Run();
private:
	// Volume renderer functionality
	void ShowImage();
	void UpdateTexture(const void* buffer, const int width, const int height);
	void Process();
	// UI
	void ControlPanel(uint32_t width, uint32_t height);
	void Scene(uint32_t width, uint32_t height);
	// SDL event related functions
	void OnResize(uint32_t width, uint32_t height);
	void HandleSDLEvent(SDL_Event* event);
	//Ray casting
	GLubyte* RayCasting();
	// SDL & window
	SDL_Window* m_window;
	SDL_GLContext m_context;
	uint32_t m_width, m_height;
	// Status
	short* pVolume;
	int vwidth;
	int vheight;
	int vvolume;
	bool m_isRunning;
	float theta = 90;
	float phi = 0;
	float tf_tail = 120;
	float tf_head = 1000;
	float trans_threshold = 0.0;
	float sampling_rate = 0.5;
	float ray_rate = 1;
	int imWidthHeight = 501;
	float kd = 0.2;
	float ka = 0.7;
	float ks = 0.14;
	int n = 100;
	// Output texture
	bool m_hasTexture;
	GLuint m_textureID;
};

