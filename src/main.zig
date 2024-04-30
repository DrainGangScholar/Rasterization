const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const CHANNELS: usize = 3;
const stdout = std.io.getStdOut();
const write = stdout.writer();
const math = std.math;

const Matrix = struct {
    buffer: [3][3]f64,
    pub fn mul(self: Matrix, other: Matrix) Matrix {
        var buffer: [3][3]f64 = undefined;
        for (0..3) |i| {
            for (0..3) |j| {
                buffer[i][j] += self.buffer[i][j] * other.buffer[j][i];
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
    for (points.items) |point| {
        index = (@abs(point.y) * WIDTH + @abs(point.x)) * CHANNELS;
        buffer[index] = 0;
        buffer[index + 1] = 128;
        buffer[index + 2] = 255;
    }
}
pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const size = CHANNELS * WIDTH * HEIGHT;
    const buffer = try allocator.alloc(u8, size);
    defer allocator.free(buffer);

    try init_buffer(buffer);
    const start = Point{ .x = 0, .y = 0, .w = 1 };
    const end = Point{ .x = WIDTH / 2 - 1, .y = HEIGHT / 2 - 1, .w = 1 };

    //    const start1 = Point{ .x = 0, .y = HEIGHT - 1, .w = 1 };
    //    const end1 = Point{ .x = WIDTH / 2 - 1, .y = HEIGHT / 2 - 1, .w = 1 };
    //
    //    const start2 = Point{ .x = 0, .y = 0, .w = 1 };
    //    const end2 = Point{ .x = 0, .y = HEIGHT - 1, .w = 1 };
    //
    //    const start3 = Point{ .x = 0, .y = HEIGHT / 2 - 1, .w = 1 };
    //    const end3 = Point{ .x = WIDTH / 2 - 1, .y = HEIGHT / 2 - 1, .w = 1 };

    //    const t = translate((WIDTH / 2 - 1), 0.0);
    //const r = rotate(45.0);
    const s = scale(-1.0, 1.0);
    const tr = s;

    try rasterize(start.mul(tr), end.mul(tr), buffer, allocator);
    //    try rasterize(start1.mul(tr), end1.mul(tr), buffer, allocator);
    //    try rasterize(start2.mul(tr), end2.mul(tr), buffer, allocator);
    //    try rasterize(start3.mul(tr), end3.mul(tr), buffer, allocator);

    try write_to_ppm(buffer);
}
