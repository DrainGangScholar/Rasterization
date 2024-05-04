const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const WIDTH: usize = 640;
const C_WIDTH: c_int = 640;
const HEIGHT: usize = 480;
const C_HEIGHT: c_int = 480;
const CHANNELS: usize = 3;
const size = CHANNELS * WIDTH * HEIGHT;
const stdout = std.io.getStdOut();
const write = stdout.writer();
const math = std.math;
const PI = std.math.pi;
const PI_OVER_180 = PI / 180.0;

const c = @cImport({
    @cInclude("SDL.h");
});

inline fn to_rad(angle: f64) f64 {
    return angle * PI_OVER_180;
}
const SDL_WINDOWPOS_UNDEFINED: c_int = @bitCast(c.SDL_WINDOWPOS_UNDEFINED_MASK);

const Matrix = struct {
    buffer: [3][3]f64,
    pub fn mul(self: Matrix, other: Matrix) Matrix {
        var buffer: [3][3]f64 = undefined;
        for (0..3) |i| {
            for (0..3) |j| {
                for (0..3) |k| {
                    buffer[i][j] += self.buffer[i][k] * other.buffer[k][j];
                }
            }
        }
        return Matrix{
            .buffer = buffer,
        };
    }
};
pub fn translate(dx: f64, dy: f64) Matrix {
    return Matrix{
        .buffer = .{
            .{ 1, 0.0, 0.0 },
            .{ 0.0, 1.0, 0.0 },
            .{ dx, dy, 1.0 },
        },
    };
}
pub fn scale(dx: f64, dy: f64) Matrix {
    return Matrix{
        .buffer = .{
            .{ dx, 0.0, 0.0 },
            .{ 0.0, dy, 0.0 },
            .{ 0.0, 0.0, 1.0 },
        },
    };
}
pub fn rotate(angle: f64) Matrix {
    const cos = math.cos(angle);
    const sin = math.sin(angle);
    return Matrix{
        .buffer = .{
            .{ cos, sin, 0.0 },
            .{ -sin, cos, 0.0 },
            .{ 0.0, 0.0, 1.0 },
        },
    };
}

const Point = struct {
    x: i64,
    y: i64,
    w: i64,
    pub fn x_usize(self: *const Point) usize {
        return @intCast(self.x);
    }
    pub fn y_usize(self: *const Point) usize {
        return @intCast(self.y);
    }
    pub fn x_float(self: *const Point) f64 {
        return @floatFromInt(self.x);
    }
    pub fn y_float(self: *const Point) f64 {
        return @floatFromInt(self.y);
    }
    pub fn w_float(self: *const Point) f64 {
        return @floatFromInt(self.w);
    }
    pub fn mul(self: *const Point, mat: Matrix) Point {
        const x = self.x_float() * mat.buffer[0][0] + self.y_float() * mat.buffer[1][0] + self.w_float() * mat.buffer[2][0];
        const y = self.x_float() * mat.buffer[0][1] + self.y_float() * mat.buffer[1][1] + self.w_float() * mat.buffer[2][1];
        const w = self.x_float() * mat.buffer[0][2] + self.y_float() * mat.buffer[1][2] + self.w_float() * mat.buffer[2][2];
        return Point{
            .x = @intFromFloat(x / w),
            .y = @intFromFloat(y / w),
            .w = 1,
        };
    }
};

pub fn BresenhamLine(start: Point, end: Point, allocator: Allocator) !ArrayList(Point) {
    var points = ArrayList(Point).init(allocator);
    errdefer points.deinit();

    const x0 = start.x;
    const y0 = start.y;
    const x1 = end.x;
    const y1 = end.y;

    const dx: i64 = @intCast(@abs(x0 - x1));
    const dy: i64 = @intCast(@abs(y0 - y1));

    var x = x0;
    var y = y0;

    const sx: i64 = if (x0 < x1) 1 else -1;
    const sy: i64 = if (y0 < y1) 1 else -1;

    var err = dx - dy;

    try points.append(Point{
        .x = x,
        .y = y,
        .w = 1,
    });

    var err2: i64 = 0;
    while (x != x1 or y != y1) {
        err2 = 2 * err;
        if (err2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (err2 < dx) {
            err += dx;
            y += sy;
        }

        try points.append(Point{
            .x = x,
            .y = y,
            .w = 1,
        });
    }
    return points;
}
pub fn write_to_ppm(buffer: []u8) !void {
    var index: usize = 0;
    try write.print("P3\n{} {}\n255\n", .{ WIDTH, HEIGHT });
    for (0..HEIGHT) |i| {
        for (0..WIDTH) |j| {
            index = (i * WIDTH + j) * CHANNELS;
            try write.print("{} {} {}\n", .{ buffer[index], buffer[index + 1], buffer[index + 2] });
        }
    }
}
pub fn init_buffer(buffer: []u8) !void {
    var index: usize = 0;
    for (0..HEIGHT) |i| {
        for (0..WIDTH) |j| {
            index = (i * WIDTH + j) * CHANNELS;
            buffer[index] = 255;
            buffer[index + 1] = 255;
            buffer[index + 2] = 255;
        }
    }
}
pub fn rasterize(start: Point, end: Point, buffer: []u8, allocator: Allocator) !void {
    const points = try BresenhamLine(start, end, allocator);
    defer points.deinit();
    var index: usize = 0;
    var x: usize = 0;
    var y: usize = 0;
    for (points.items) |point| {
        x = if (point.x < 0) 0 else if (point.x >= WIDTH) WIDTH - 1 else @intCast(point.x);
        y = if (point.y < 0) 0 else if (point.y >= HEIGHT) HEIGHT - 1 else @intCast(point.y);
        index = (y * WIDTH + x) * CHANNELS;
        buffer[index] = 0;
        buffer[index + 1] = 128;
        buffer[index + 2] = 255;
    }
}
pub fn main() !void {
    if (c.SDL_Init(c.SDL_INIT_VIDEO) != 0) {
        c.SDL_Log("Unable to initialize SDL: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    }
    defer c.SDL_Quit();
    const window = c.SDL_CreateWindow("Paint", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, C_WIDTH, C_HEIGHT, c.SDL_WINDOW_OPENGL | c.SDL_WINDOW_RESIZABLE) orelse
        {
        c.SDL_Log("Unable to create window: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    };

    const renderer = c.SDL_CreateRenderer(window, -1, 0) orelse {
        c.SDL_Log("Unable to create renderer: %s", c.SDL_GetError());
        return error.SDLInitializationFailed;
    };
    defer c.SDL_DestroyRenderer(renderer);

    defer c.SDL_DestroyWindow(window);
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const buffer = try allocator.alloc(u8, size);
    defer allocator.free(buffer);

    try init_buffer(buffer);
    const side_length = 60;

    const center_x: i64 = @intCast(WIDTH / 2);
    const center_y: i64 = @intCast(HEIGHT / 2);

    const half_height: i64 = @intFromFloat(side_length * math.sqrt(3.0) / 2.0);
    const half_width: i64 = @intFromFloat(side_length / 2.0);

    const start = Point{ .x = center_x, .y = center_y + half_height, .w = 1 };
    const end = Point{ .x = center_x, .y = center_y - half_height, .w = 1 };

    const start1 = Point{ .x = center_x - half_width, .y = center_y + half_height, .w = 1 };
    const end1 = Point{ .x = center_x + half_width, .y = center_y + half_height, .w = 1 };

    const start2 = Point{ .x = center_x - half_width, .y = center_y + half_height, .w = 1 };
    const end2 = Point{ .x = center_x, .y = center_y - half_height, .w = 1 };

    const start3 = Point{ .x = center_x + half_width, .y = center_y + half_height, .w = 1 };
    const end3 = Point{ .x = center_x, .y = center_y - half_height, .w = 1 };

    const translate_to_origin = translate(-center_x, -center_y);
    const translate_back = translate(center_x, center_y);

    const angle = to_rad(0.0);
    const r = rotate(angle);
    const s = scale(1.0, -1.0);

    const trt = translate_to_origin.mul(r.mul(s.mul(translate_back)));

    try rasterize(start.mul(trt), end.mul(trt), buffer, allocator);
    try rasterize(start1.mul(trt), end1.mul(trt), buffer, allocator);
    try rasterize(start2.mul(trt), end2.mul(trt), buffer, allocator);
    try rasterize(start3.mul(trt), end3.mul(trt), buffer, allocator);

    var quit: bool = false;
    while (!quit) {
        var event: c.SDL_Event = undefined;
        while (c.SDL_PollEvent(&event) != 0) {
            switch (event.type) {
                c.SDL_KEYDOWN => {
                    const keycode = event.key.keysym.sym;
                    if (keycode == c.SDLK_ESCAPE) {
                        quit = true;
                    }
                },
                c.SDL_MOUSEBUTTONDOWN => {
                    if (event.button.button == c.SDL_BUTTON_LEFT) {
                        try write.print("MOUSE CLICK\n", .{});
                    }
                },
                else => {},
            }
        }

        c.SDL_Delay(10);
    }
    //    try write_to_ppm(buffer);
}
