#pragma once
#include <string>
#include <GL/glew.h>
class TextureSky {
private:
	unsigned int m_RendererID;
	std::string m_FilePath;
	unsigned char* m_LocalBuffer;
	int m_Width, m_Height, m_BPP;

public:
	TextureSky(const std::string& path);
	~TextureSky();

	void Bind(unsigned int slot) const;
	void Unbind() const;

	inline int GetWidth() { return m_Width; }
	inline int GetHeight() { return m_Height; }

};