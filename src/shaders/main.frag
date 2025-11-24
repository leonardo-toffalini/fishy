#version 330

in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;

uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform float light_intensity;
uniform vec3 colormap[256];
uniform int colormap_size;
uniform int reversed;

out vec4 finalColor;

void main() {
  vec3 lightDir = normalize(vec3(1.0, 1.0, 1.0));
  vec3 lightColor = vec3(1.0, 1.0, 1.0);
  vec3 ambient = vec3(0.3, 0.3, 0.3);

  float diff = max(dot(normalize(fragNormal), lightDir), 0.0);
  vec3 diffuse = diff * lightColor;

  vec3 lighting = ambient + diffuse;

  vec4 texelColor = texture(texture0, fragTexCoord);
  float val;
  if (reversed == 0)
    val = clamp(texelColor.r, 0.0, 1.0);
  else
    val = 1 - clamp(texelColor.r, 0.0, 1.0);

  int last_color_idx = colormap_size - 1;
  float segments = float(last_color_idx);
  float scaled = val * segments;
  int index = int(floor(scaled));
  float t = fract(scaled);

  vec3 baseColor = colormap[index];
  vec3 nextColor = colormap[min(index + 1, last_color_idx)];
  vec3 color = mix(baseColor, nextColor, t);

  if (light_intensity > 0.5)
    finalColor = vec4(color * lighting, texelColor.a);
  else
    finalColor = vec4(color, texelColor.a);
}


