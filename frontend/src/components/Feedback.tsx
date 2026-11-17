import type { ReactNode } from "react";

export function Alert({
  kind,
  children,
}: {
  kind: "info" | "warning" | "error";
  children: ReactNode;
}) {
  return (
    <div className={`alert alert-${kind}`} role={kind === "error" ? "alert" : "status"}>
      {children}
    </div>
  );
}

export function Spinner({ children }: { children: ReactNode }) {
  return (
    <div className="spinner" role="status">
      <span className="spinner-ring" aria-hidden="true" />
      <em>{children}</em>
    </div>
  );
}

export function Expander({ title, children }: { title: string; children: ReactNode }) {
  return (
    <details className="expander">
      <summary>{title}</summary>
      <div className="expander-body">{children}</div>
    </details>
  );
}
