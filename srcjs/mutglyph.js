import { embed } from "@genome-spy/core/minimal";

HTMLWidgets.widget({
  name: "mutglyph",
  type: "output",

  factory(el) {
    let api;
    let renderId = 0;

    el.style.position = "relative";
    el.style.overflow = "hidden";

    return {
      async renderValue(x) {
        const currentRenderId = ++renderId;

        api?.finalize();
        api = undefined;

        const container = document.createElement("div");
        container.style.width = "100%";
        container.style.height = "100%";
        el.replaceChildren(container);

        const nextApi = await embed(container, x.spec);

        if (currentRenderId === renderId) {
          api = nextApi;
        } else {
          nextApi.finalize();
          container.remove();
        }
      },

      resize() {
        // GenomeSpy observes its container with ResizeObserver.
      },
    };
  },
});
