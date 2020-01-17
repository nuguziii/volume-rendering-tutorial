#include "Editor.h"
#include "Vector.h"
#include <iostream>
#include <regex>

#include <fstream>
#include <math.h>
#include <intrin.h>
#define _CRT_SECURE_NO_WARNINGS
#include <climits>

void TransferFunc(const float& intensity, const float& tf_head, const float& tf_tail, float& val);
void Interpolation(const float& xx, const float& yy, const float& zz, short& a, short& b, short& c, short& d, short& e, short& f, short& g, short& h, float& intensity, Vec& normVec);
float getNormVector(float a, float b, float interval);
void GetPVolumeVal(const short* pVolume, const int& width, const int& height, const int& volume, const Vec& scaled_target, short& a, short& b, short& c, short& d, short& e, short& f, short& g, short& h);

int16_t swap_int16(int16_t val){
	return (val << 8) | ((val >> 8) & 0xFF);
}

Editor::Editor(uint32_t width, uint32_t height)
	:m_window(nullptr),
	m_context(nullptr),
	m_width(width),
	m_height(height),
	m_isRunning(true),
	m_hasTexture(false) {
}

Editor::~Editor() {
	ImGui_ImplSdlGL3_Shutdown();
	SDL_GL_DeleteContext(m_context);
	SDL_DestroyWindow(m_window);
	SDL_Quit();
}

bool Editor::Initialize() {
	// Setup SDL
	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
		std::cout << ("Error: %s\n", SDL_GetError()) << std::endl;
		return false;
	}

	// Setup window
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
	SDL_DisplayMode current;
	SDL_GetCurrentDisplayMode(0, &current);
	m_window = SDL_CreateWindow("Volume Renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, m_width, m_height, SDL_WINDOW_OPENGL);
	SDL_GLContext glcontext = SDL_GL_CreateContext(m_window);
	glewInit();

	// Setup ImGui binding
	ImGui_ImplSdlGL3_Init(m_window);
	Process();
	return true; // Return initialization result
}

void Editor::Run() {
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
	while (m_isRunning) {
		// Handle SDL events
		SDL_Event event;
		if (SDL_PollEvent(&event)) {
			ImGui_ImplSdlGL3_ProcessEvent(&event);
			HandleSDLEvent(&event);
		}
		ImGui_ImplSdlGL3_NewFrame(m_window);
		// Editor
		{
			ControlPanel(m_width - 720, 720);
			Scene(720, 720);
			// Code sample of ImGui (Remove comment when you want to see it)
			//ImGui::ShowTestWindow();
		}
		// Rendering
		glViewport(0, 0, (int)ImGui::GetIO().DisplaySize.x, (int)ImGui::GetIO().DisplaySize.y);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);

		ImGui::Render();
		SDL_GL_SwapWindow(m_window);
	}
	delete[] pVolume;
}

void Editor::UpdateTexture(const void * buffer, int width, int height) {
	if (!m_hasTexture) {
		auto err = glGetError();
		glGenTextures(1, &m_textureID);
		if (err != GL_NO_ERROR) {
			throw std::runtime_error("Not able to create texture from buffer" + std::to_string(glGetError()));
		}
		else {
			m_hasTexture = true;
		}
	}
	glBindTexture(GL_TEXTURE_2D, m_textureID);
	// set texture sampling methods
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
	glBindTexture(GL_TEXTURE_2D, 0);
}

float getNormVector(float a, float b, float interval) {
	float left = b * (interval - 0.01) + a * (1 - interval + 0.01);
	float right = b * (interval + 0.01) + a * (1 - interval - 0.01);
	return right - left;
}

void Interpolation(const float& xx, const float& yy, const float& zz, short& a, short& b, short& c, short& d, short& e, short& f, short& g, short& h, float & intensity, Vec &normVec) {
	//Trilinear interpolation
	float x = xx - floor(xx);
	float y = yy - floor(yy);
	float z = zz - floor(zz);

	float ae = e * x + (1 - x) * a;
	float cg = g * x + (1 - x) * c;
	float aecg = ae * z + (1 - z) * cg;
	float bf = f * x + (1 - x) * b;
	float dh = h * x + (1 - x) * d;
	float bfdh = bf * z + (1 - z) * dh;
	normVec.y = getNormVector(aecg, bfdh, y);
	intensity = bfdh * y + (1 - y) * aecg;

	float ab = b * y + a * (1 - y);
	float ef = f * y + e * (1 - y);
	float abef = ef * x + ab * (1 - x);
	float cd = d * y + c * (1 - y);
	float gh = h * y + g * (1 - y);
	float cdgh = gh * x + cd * (1 - x);
	normVec.z = getNormVector(cdgh, abef, z);

	float abcd = ab * z + cd * (1 - z);
	float efgh = ef * z + gh * (1 - z);
	normVec.x = getNormVector(abcd, efgh, x);

	normVec = normVec.normalised();
}

void TransferFunc(const float &intensity, const float &tf_head, const float &tf_tail, float &val) {
	if (intensity <= tf_tail) val = 0;
	else if (intensity >= tf_head) val = 255;
	else {
		float cof = (255 / (tf_head - tf_tail));
		val = (cof * (intensity - (tf_tail)));
	}
}

void GetPVolumeVal(const short* pVolume, const int &width, const int &height, const int &volume, const Vec &scaled_target, short& a, short& b, short& c, short& d, short& e, short& f, short& g, short& h) {
	int idx;
	Vec coord;
	coord = Vec(floor(scaled_target.x), floor(scaled_target.y), ceil(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + width * (height - 1 - coord.z) + coord.y);
	a = pVolume[idx];
	coord = Vec(floor(scaled_target.x), ceil(scaled_target.y), ceil(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	b = pVolume[idx];
	coord = Vec(floor(scaled_target.x), floor(scaled_target.y), floor(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	c = pVolume[idx];
	coord = Vec(floor(scaled_target.x), ceil(scaled_target.y), floor(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	d = pVolume[idx];
	coord = Vec(ceil(scaled_target.x), floor(scaled_target.y), ceil(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	e = pVolume[idx];
	coord = Vec(ceil(scaled_target.x), ceil(scaled_target.y), ceil(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	f = pVolume[idx];
	coord = Vec(ceil(scaled_target.x), floor(scaled_target.y), floor(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	g = pVolume[idx];
	coord = Vec(ceil(scaled_target.x), ceil(scaled_target.y), floor(scaled_target.z));
	idx = static_cast<int>(width * height * (volume - 1 - coord.x) + (width) * (height - 1 - coord.z) + coord.y);
	h = pVolume[idx];
}

GLubyte* Editor::RayCasting() {
	float width = (vwidth - 1) * 0.66;
	float height = (vheight - 1) * 0.66;
	float volume = (vvolume - 1) * 5;

	GLubyte* imVals = new GLubyte[imWidthHeight * imWidthHeight * 4];
	const float sinTheta = abs(sin((theta) * M_PI / 180)) < 1e-5 ? 0 : sin((theta) * M_PI / 180);
	const float cosTheta = abs(cos((theta) * M_PI / 180)) < 1e-5 ? 0 : cos((theta) * M_PI / 180);
	const float sinPhi = abs(sin(phi * M_PI / 180)) < 1e-5 ? 0 : sin(phi * M_PI / 180);
	const float cosPhi = abs(cos(phi * M_PI / 180)) < 1e-5 ? 0 : cos(phi * M_PI / 180);
	float radius = sqrt(pow(volume, 2) + pow(width, 2) + pow(height, 2))/2.0;

	const float LsinTheta = abs(sin((60 + theta) * M_PI / 180)) < 1e-5 ? 0 : sin((60 + theta) * M_PI / 180);
	const float LcosTheta = abs(cos((60 + theta) * M_PI / 180)) < 1e-5 ? 0 : cos((60 + theta) * M_PI / 180);
	const float LsinPhi = abs(sin((90 + phi) * M_PI / 180)) < 1e-5 ? 0 : sin((90 + phi) * M_PI / 180);
	const float LcosPhi = abs(cos((90 + phi) * M_PI / 180)) < 1e-5 ? 0 : cos((90 + phi) * M_PI / 180);

	Vec center_coordinate = Vec(radius * sinTheta * cosPhi, radius * cosTheta, radius * sinTheta * sinPhi);
	Vec light_coordinate = Vec(radius * LsinTheta * LcosPhi, radius * LcosTheta, radius * LsinTheta * LsinPhi);
	Vec L = light_coordinate.normalised() * -1.0;
	Vec V = center_coordinate.normalised() * -1.0;
	Vec leftVec = Vec(-1* V.y, V.x, 0);
	leftVec = leftVec.normalised();
	Vec upVec = Vec(-1*leftVec.y * V.z, leftVec.x * V.z, leftVec.y * V.x - leftVec.x * V.y);
	upVec = upVec.normalised();

	Vec start_coordinate = center_coordinate+(upVec + leftVec)*static_cast<int>(imWidthHeight /2)* ray_rate;

	for (int i = 0;i < imWidthHeight;i++) {
		for (int j = 0;j < imWidthHeight;j++) {
			Vec cur_coordinate = start_coordinate + (leftVec * -1 * j + upVec * -1 * i)* ray_rate;
			float fVal = 0.0;
			float tVal = 1.0;
			//ray casting
			for (float length = 1;length < radius*2;length += sampling_rate) {
				Vec cur_ray_coordinate = cur_coordinate + V * length;
				Vec bottom_range = Vec(-1 * (volume / 2), -1 * (width / 2), -1 * (height / 2));
				Vec top_range = Vec((volume / 2), (width / 2), (height / 2));

				Vec center_ray_coordinate = cur_coordinate + V * radius;
				if (!(center_ray_coordinate.x >= bottom_range.x && center_ray_coordinate.x <= top_range.x &&
					center_ray_coordinate.y >= bottom_range.y && center_ray_coordinate.y <= top_range.y &&
					center_ray_coordinate.z >= bottom_range.z && center_ray_coordinate.z <= top_range.z)) {
					break;
				}

				if ((cur_ray_coordinate.x>=bottom_range.x) && (cur_ray_coordinate.x<= top_range.x) &&
					(cur_ray_coordinate.y >= bottom_range.y) && (cur_ray_coordinate.y <= top_range.y) &&
					(cur_ray_coordinate.z >= bottom_range.z) && (cur_ray_coordinate.z <= top_range.z)){

					Vec N = Vec();
					Vec H = (V + L).normalised();

					float ambient = ka;

					//interpolation
					short a, b, c, d, e, f, g, h;
					float intensity;
					Vec scaled_target = Vec((volume/2.0 + cur_ray_coordinate.x) / 5.0, ((width/2.0) + cur_ray_coordinate.y) / 0.66, (height/2.0 + cur_ray_coordinate.z)/ 0.66);
					GetPVolumeVal(pVolume, vwidth, vheight, vvolume, scaled_target, a, b, c, d, e, f, g, h);
					Interpolation(scaled_target.x, scaled_target.y, scaled_target.z, a, b, c, d, e, f, g, h, intensity, N);

					//transfer function
					float val = 0.0;
					TransferFunc(intensity, tf_head, tf_tail, val);

					float diffuseLight = L ^ N;
					diffuseLight = diffuseLight > 0 ? diffuseLight : 0;
					float diffuse = kd * diffuseLight;

					float HN = H ^ V;
					HN = HN > 0 ? HN : 0;
					float specularLight = pow(HN, n);

					if (diffuseLight < 0)
						specularLight = 0;

					float specular = ks * specularLight;

					//front-to-back compositing update
					float shading = (val / 255) * (ambient + diffuse + specular);
					fVal = fVal + (val / 255) * shading * tVal;
					tVal = (1 - (val / 255)) * tVal;
					if (tVal <= trans_threshold) break;
				}
			}
			//phong shader
			//Vec light = Vec(sinTheta * cosPhi, cosTheta, sinTheta * (abs(sin((phi+90) * M_PI / 180)) < 1e-5 ? 0 : sin(phi * M_PI / 180))).normalised();
			//fVal = light * kd * fVal + (light + dirVec) * dirVec * ks;
			imVals[4 * (imWidthHeight * i + j)] = static_cast<int>(255 * fVal);
			imVals[4 * (imWidthHeight * i + j)+1] = static_cast<int>(255 * fVal);
			imVals[4 * (imWidthHeight * i + j)+2] = static_cast<int>(255 * fVal);
			imVals[4 * (imWidthHeight * i + j) + 3] = static_cast<int>(255 * (1-tVal));
		}
	}
	return imVals;
}

void Editor::ShowImage() {
	std::cout << "ray casting...";
	GLubyte* imVals = RayCasting();

	std::cout << "done!" << std::endl;
	UpdateTexture(imVals, imWidthHeight, imWidthHeight);

	delete[] imVals;
}

void Editor::Process() {
	/* TODO : Process volume data & pass raw buffer to UpdateTexture method*/
	vwidth = 512;
	vheight = 512;
	vvolume = 56;

	const int volSize = vwidth * vheight * vvolume;

	//load file data
	const char* fileName = "../asset/data/volume1_512x512x56-short-bigendian.raw";
	pVolume = new short[volSize];

	FILE* pFile;
	errno_t	err = fopen_s(&pFile, fileName, "rb");

	if (err!=0 || !pFile) {
		std::cerr << "Error opening file" << std::endl;
		std::exit(EXIT_FAILURE);
	}

	fread(pVolume, sizeof(short), volSize, pFile);
	fclose(pFile);

	for (int i = 0;i < volSize;i++) {
		pVolume[i] = swap_int16(pVolume[i]);
	}
	std::cout << "load file...done!" << std::endl;

	ShowImage();
}

void Editor::ControlPanel(uint32_t width, uint32_t height) {
	// Control Panel Window
	ImGui::SetNextWindowSize(ImVec2((float)width, (float)height));
	ImGui::SetNextWindowPos(ImVec2(0, 0));
	ImGui::Begin("Control Panel", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);

	/* TODO : Write UI Functions */
	ImGui::Text("direction");
	ImGui::SliderFloat("up-down", &phi, 0.0f, 359.99f);
	ImGui::SliderFloat("left-right", &theta, 0.0f, 359.99f);
	ImGui::Text("");
	ImGui::Text("transfer function");
	ImGui::SliderFloat("bottom", &tf_tail, -2000.0f, 2000.0f);
	ImGui::SliderFloat("top", &tf_head, -2000.0f, 2000.0f);
	ImGui::Text("");
	ImGui::Text("transparency (Tn)");
	ImGui::SliderFloat("", &trans_threshold, 0.0f, 1.0f);
	ImGui::Text("");
	ImGui::Text("sampling/ray interval");
	ImGui::SliderFloat("sampling", &sampling_rate, 0.03, 2);
	ImGui::SliderFloat("ray", &ray_rate, 0.1, 1);
	ImGui::Text("");
	ImGui::Text("phong shading params");
	ImGui::SliderFloat("ambient", &ka, 0, 1);
	ImGui::SliderFloat("diffuse", &kd, 0, 1);
	ImGui::SliderFloat("specular", &ks, 0, 1);
	ImGui::SliderInt("shininess", &n, 1, 256);

	ImGui::End();
}

void Editor::Scene(uint32_t width, uint32_t height) {
	// Scene Window
	ImGui::SetNextWindowSize(ImVec2((float)width, (float)height));
	ImGui::SetNextWindowPos(ImVec2((float)(m_width - width), 0.f));
	ImGui::Begin("Scene", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
	// Draw texture if there is one
	if (m_hasTexture) {
		ImGui::Image(ImTextureID(m_textureID), ImGui::GetContentRegionAvail());
	}
	ImGui::End();
}

void Editor::OnResize(uint32_t width, uint32_t height) {
	m_width = width;
	m_height = height;
}

void Editor::HandleSDLEvent(SDL_Event * event) {
	// SDL_Event wiki : https://wiki.libsdl.org/SDL_Event
	static bool mouseIsDown = false;
	static bool isDragging = false;
	short moved = event->motion.xrel;
	int degreeStep = 5;
	//short interval = abs(moved) / 5;
	short interval = 1;
	switch (event->type) {
	case SDL_QUIT:
		m_isRunning = false;
		break;
	case SDL_KEYDOWN:
		break;
	case SDL_MOUSEWHEEL:
		break;
	case SDL_MOUSEMOTION:
		if (mouseIsDown) {
			ShowImage();
		}
		break;
	case SDL_MOUSEBUTTONDOWN:
		mouseIsDown = true;
		ShowImage();
		break;
	case SDL_MOUSEBUTTONUP:
		mouseIsDown = false;
		break;
	case SDL_WINDOWEVENT:
		switch (event->window.event) {
		case SDL_WINDOWEVENT_RESIZED:
			OnResize(event->window.data1, event->window.data2);
			break;
		default:
			break;
		}
	default:
		break;
	}
}
