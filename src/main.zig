const std = @import("std");
const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const Point = struct {
    x: i32,
    y: i32,
};
pub fn BresenhamLine(start: Point, end: Point, allocator: Allocator) !ArrayList(Point) {
    var points = ArrayList(Point).init(allocator);
    errdefer points.deinit();
    const dx: i32 = @intCast(@abs(start.x - end.x));
    const dy: i32 = @intCast(@abs(start.y - end.y));

    const start_greater: bool = start.x > end.x;
    var x = if (start_greater) end.x else start.x;
    var y = if (start_greater) end.y else start.y;
    const xend = if (start_greater) start.x else end.x;

    var d = 2 * dy - dx;
    const incr1: i32 = 2 * dy;
    const incr2: i32 = 2 * (dy - dx);

    try points.append(Point{
        .x = x,
        .y = y,
    });

    while (x < xend) {
        x += 1;
        if (d < 0) {
            d += incr1;
        } else {
            y += 1;
            d += incr2;
        }
        try points.append(Point{
            .x = x,
            .y = y,
        });
    }
    return points;
}
pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();
    const start = Point{ .x = 0, .y = 0 };
    const end = Point{ .x = 5, .y = 3 };
    const points = try BresenhamLine(start, end, allocator);
    defer points.deinit();
    for (points.items) |point| {
        std.debug.print("{}\n", .{point});
    }
}
