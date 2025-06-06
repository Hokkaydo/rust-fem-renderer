use glyphon::{Attrs, Buffer, Color, Family, FontSystem, Metrics, Shaping, TextArea, TextBounds};
use winit::window::Window;


pub struct Texts {
    texts: Vec<Text>,
}

pub struct Text {
    id: String,
    text_buffer: Buffer,
    x: f32,
    y: f32,
    color: Color,
    visible: bool,
}

impl Texts {
    pub fn new() -> Self {
        Self { texts: Vec::new() }
    }

    pub fn add_text(&mut self, text: Text) {
        self.texts.push(text);
    }

    pub fn render(&self, window: &Window) -> Vec<TextArea> {
        self.texts.iter().filter(|text| text.visible).map(|text| text.render(window)).collect()
    }

    pub fn get_text(&mut self, id: &str) -> Option<&mut Text> {
        self.texts.iter_mut().find(|text| text.id == id)
    }

    pub fn fonts_resize(&mut self, font_system: &mut FontSystem, size: f32) {
        for text in &mut self.texts {
            text.set_font_size(font_system, size);
        }
    }

}

impl Text {

    pub fn new(
        id: &str,
        content: &str,
        font_system: &mut FontSystem,
        window: &Window,
        x: f32,
        y: f32,
        color: Color,
        visible: bool,
    ) -> Self {
        let mut text_buffer = Buffer::new(font_system, Metrics::new(24.0, 24.0));

        let physical_width = (window.inner_size().width as f64 * window.scale_factor()) as f32;
        let physical_height = (window.inner_size().height as f64 * window.scale_factor()) as f32;

        text_buffer.set_size(font_system, Some(physical_width), Some(physical_height));
        text_buffer.set_text(
            font_system,
            content,
            Attrs::new().family(Family::SansSerif),
            Shaping::Advanced,
        );
        text_buffer.shape_until_scroll(font_system, false);

        Self { id: id.to_string(), text_buffer, x, y, color, visible }
    }

    pub fn set_text(&mut self, content: &str, font_system: &mut FontSystem) {
        self.text_buffer.set_text(
            font_system,
            content,
            Attrs::new().family(Family::SansSerif),
            Shaping::Advanced,
        );
        self.text_buffer.shape_until_scroll(font_system, false);
    }

    pub fn set_font_size(&mut self, font_system: &mut FontSystem, size: f32) {
        self.text_buffer.set_metrics(font_system, Metrics::new(size, size));
        self.text_buffer.shape_until_scroll(font_system, false);
    }

    pub fn render(&self, window: &Window) -> TextArea {
        let left =
            self.x * (window.inner_size().width as f32 * window.scale_factor() as f32) as f32;
        let top =
            self.y * (window.inner_size().height as f32 * window.scale_factor() as f32) as f32;

        TextArea {
            buffer: &self.text_buffer,
            left,
            top,
            scale: 1.0,
            bounds: TextBounds::default(),
            default_color: self.color,
            custom_glyphs: &[],
        }
    }

    pub fn toggle_visibility(&mut self) {
        self.visible = !self.visible;
    }
}
